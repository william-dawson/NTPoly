// bench_cusparse.cu  — v2
// Standalone GPU SpGEMM benchmark using cuSPARSE (CUDA 13.1 API).
// Tries multiple SpGEMM algorithms if the first fails.
//
// Build:  nvcc -O2 -std=c++17 -o bench_cusparse bench_cusparse.cu -lcusparse
// Run:    ./bench_cusparse matrix.bin [--rounds 10] [--warmup 3] [--threshold 1e-4] [--alg N]
//
// NTPoly binary format:
//   header: int32 rows, int32 cols, int32 complex_flag, int64 total_nnz  (20 bytes)
//   body:   N × {int32 col, int32 row, double64 val}  (16 bytes each)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

#include <cusparse.h>
#include <cuda_runtime.h>

#define CUDA_CHECK(call) do {                                          \
    cudaError_t err = (call);                                           \
    if (err != cudaSuccess) {                                           \
        fprintf(stderr, "CUDA error at %s:%d: %s\n",                    \
                __FILE__, __LINE__, cudaGetErrorString(err));           \
        exit(1);                                                        \
    }                                                                   \
} while(0)

// Read NTPoly binary matrix into COO format
struct COOMatrix {
    int rows, cols, nnz;
    std::vector<int>    coo_row, coo_col;
    std::vector<double> coo_val;
};

static COOMatrix read_ntpoly_binary(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); exit(1); }
    int header[3];
    long long total_nnz;
    if (fread(header, sizeof(int), 3, fp) != 3) { fprintf(stderr, "read err\n"); exit(1); }
    if (fread(&total_nnz, sizeof(long long), 1, fp) != 1) { fprintf(stderr, "read err\n"); exit(1); }
    COOMatrix mat;
    mat.rows = header[0]; mat.cols = header[1]; mat.nnz = (int)total_nnz;
    fprintf(stderr, "Matrix: %d x %d, nnz = %d\n", mat.rows, mat.cols, mat.nnz);
    mat.coo_row.resize(mat.nnz); mat.coo_col.resize(mat.nnz); mat.coo_val.resize(mat.nnz);
    for (int i = 0; i < mat.nnz; i++) {
        int col, row; double val;
        if (fread(&col, sizeof(int), 1, fp) != 1) { fprintf(stderr, "read err\n"); exit(1); }
        if (fread(&row, sizeof(int), 1, fp) != 1) { fprintf(stderr, "read err\n"); exit(1); }
        if (fread(&val, sizeof(double), 1, fp) != 1) { fprintf(stderr, "read err\n"); exit(1); }
        mat.coo_col[i] = col - 1;  // 1-based to 0-based
        mat.coo_row[i] = row - 1;
        mat.coo_val[i] = val;
    }
    fclose(fp);
    return mat;
}

// Convert COO to CSR (sorted by row, then by column within each row)
struct CSRMatrix {
    int rows, cols, nnz;
    std::vector<int>    csr_row_ptr, csr_col_ind;
    std::vector<double> csr_val;
};

static CSRMatrix coo_to_csr(const COOMatrix &coo) {
    CSRMatrix csr;
    csr.rows = coo.rows; csr.cols = coo.cols; csr.nnz = coo.nnz;
    csr.csr_row_ptr.assign(coo.rows + 1, 0);
    csr.csr_col_ind.resize(coo.nnz); csr.csr_val.resize(coo.nnz);
    for (int i = 0; i < coo.nnz; i++) csr.csr_row_ptr[coo.coo_row[i] + 1]++;
    for (int i = 0; i < coo.rows; i++) csr.csr_row_ptr[i + 1] += csr.csr_row_ptr[i];
    std::vector<int> offset(coo.rows, 0);
    for (int i = 0; i < coo.nnz; i++) {
        int r = coo.coo_row[i];
        int pos = csr.csr_row_ptr[r] + offset[r]++;
        csr.csr_col_ind[pos] = coo.coo_col[i];
        csr.csr_val[pos]     = coo.coo_val[i];
    }
    // Sort each row by column index
    for (int r = 0; r < coo.rows; r++) {
        int s = csr.csr_row_ptr[r], e = csr.csr_row_ptr[r + 1];
        std::vector<int> idx(e - s);
        for (int i = 0; i < e - s; i++) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return csr.csr_col_ind[s + a] < csr.csr_col_ind[s + b]; });
        std::vector<int> tc(e - s); std::vector<double> tv(e - s);
        for (int i = 0; i < e - s; i++) { tc[i] = csr.csr_col_ind[s+idx[i]]; tv[i] = csr.csr_val[s+idx[i]]; }
        for (int i = 0; i < e - s; i++) { csr.csr_col_ind[s+i] = tc[i]; csr.csr_val[s+i] = tv[i]; }
    }
    return csr;
}

struct BenchResult {
    double median_ms, min_ms, max_ms;
    int    result_nnz, nnz_above_threshold, nnz_per_row;
    double trace, grand_sum;
    int    alg_used;
    bool   success;
};

// Try one SpGEMM algorithm. Returns true on success.
static bool try_spgemm(cusparseHandle_t handle,
                       cusparseSpMatDescr_t desc_A,
                       int A_rows, int A_cols,
                       int *d_A_row_ptr,  // not used directly
                       int *d_C_row_ptr,
                       double alpha_val, double beta_val,
                       cusparseSpGEMMAlg_t alg,
                       int warmup, int rounds,
                       int &last_C_nnz,
                       int *&d_C_col_ind, double *&d_C_val,
                       std::vector<double> &times) {
    d_C_col_ind = nullptr;
    d_C_val = nullptr;
    times.clear();
    last_C_nnz = 0;

    for (int round = 0; round < warmup + rounds; round++) {
        if (d_C_col_ind) { cudaFree(d_C_col_ind); d_C_col_ind = nullptr; }
        if (d_C_val)     { cudaFree(d_C_val);     d_C_val = nullptr;     }

        cusparseSpMatDescr_t desc_C;
        cusparseCreateCsr(&desc_C, A_rows, A_cols, 0,
            d_C_row_ptr, nullptr, nullptr,
            CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);

        cusparseSpGEMMDescr_t spgemmDescr;
        cusparseSpGEMM_createDescr(&spgemmDescr);

        double alpha = alpha_val, beta = beta_val;
        cudaEvent_t start, stop;
        bool timing = (round >= warmup);
        if (timing) { cudaEventCreate(&start); cudaEventCreate(&stop); cudaEventRecord(start); }

        // Step 1: work estimation
        size_t bufSize1 = 0;
        cusparseStatus_t s1a = cusparseSpGEMM_workEstimation(
            handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, desc_A, desc_A, &beta, desc_C,
            CUDA_R_64F, alg, spgemmDescr, &bufSize1, nullptr);
        if (s1a != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "  workEstimation(query) alg=%d: %s\n", (int)alg, cusparseGetErrorString(s1a));
            cusparseSpGEMM_destroyDescr(spgemmDescr);
            cusparseDestroySpMat(desc_C);
            return false;
        }
        void *buf1 = nullptr;
        if (bufSize1 > 0) CUDA_CHECK(cudaMalloc(&buf1, bufSize1));
        cusparseStatus_t s1b = cusparseSpGEMM_workEstimation(
            handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, desc_A, desc_A, &beta, desc_C,
            CUDA_R_64F, alg, spgemmDescr, &bufSize1, buf1);
        if (s1b != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "  workEstimation alg=%d: %s\n", (int)alg, cusparseGetErrorString(s1b));
            if (buf1) cudaFree(buf1);
            cusparseSpGEMM_destroyDescr(spgemmDescr);
            cusparseDestroySpMat(desc_C);
            return false;
        }

        // Step 2: compute
        size_t bufSize2 = 0;
        cusparseStatus_t s2a = cusparseSpGEMM_compute(
            handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, desc_A, desc_A, &beta, desc_C,
            CUDA_R_64F, alg, spgemmDescr, &bufSize2, nullptr);
        if (s2a != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "  compute(query) alg=%d: %s\n", (int)alg, cusparseGetErrorString(s2a));
            if (buf1) cudaFree(buf1);
            cusparseSpGEMM_destroyDescr(spgemmDescr);
            cusparseDestroySpMat(desc_C);
            return false;
        }
        void *buf2 = nullptr;
        if (bufSize2 > 0) CUDA_CHECK(cudaMalloc(&buf2, bufSize2));
        cusparseStatus_t s2b = cusparseSpGEMM_compute(
            handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, desc_A, desc_A, &beta, desc_C,
            CUDA_R_64F, alg, spgemmDescr, &bufSize2, buf2);
        if (s2b != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "  compute alg=%d: %s\n", (int)alg, cusparseGetErrorString(s2b));
            if (buf1) cudaFree(buf1);
            if (buf2) cudaFree(buf2);
            cusparseSpGEMM_destroyDescr(spgemmDescr);
            cusparseDestroySpMat(desc_C);
            return false;
        }

        int64_t C_rows, C_cols, C_nnz;
        cusparseSpMatGetSize(desc_C, &C_rows, &C_cols, &C_nnz);
        last_C_nnz = (int)C_nnz;

        CUDA_CHECK(cudaMalloc(&d_C_col_ind, sizeof(int) * C_nnz));
        CUDA_CHECK(cudaMalloc(&d_C_val, sizeof(double) * C_nnz));
        cusparseCsrSetPointers(desc_C, d_C_row_ptr, d_C_col_ind, d_C_val);

        // Step 3: copy
        cusparseStatus_t s3 = cusparseSpGEMM_copy(
            handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, desc_A, desc_A, &beta, desc_C,
            CUDA_R_64F, alg, spgemmDescr);
        if (s3 != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "  copy alg=%d: %s\n", (int)alg, cusparseGetErrorString(s3));
            if (buf1) cudaFree(buf1);
            if (buf2) cudaFree(buf2);
            cusparseSpGEMM_destroyDescr(spgemmDescr);
            cusparseDestroySpMat(desc_C);
            return false;
        }

        cudaDeviceSynchronize();

        if (timing) {
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            float ms = 0;
            cudaEventElapsedTime(&ms, start, stop);
            times.push_back(ms);
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
        }

        if (buf1) cudaFree(buf1);
        if (buf2) cudaFree(buf2);
        cusparseSpGEMM_destroyDescr(spgemmDescr);

        // Keep last desc_C for result extraction on final round
        if (round < warmup + rounds - 1) {
            cusparseDestroySpMat(desc_C);
        }
        // desc_C from the last round is leaked here but we extract data below;
        // actually let's just extract from d_C_* arrays directly.
        cusparseDestroySpMat(desc_C);
    }
    return true;
}

static BenchResult benchmark_spgemm(const CSRMatrix &A,
                                     int warmup, int rounds,
                                     double threshold,
                                     int alg_hint) {
    BenchResult res = {};
    res.success = false;
    res.median_ms = -1;

    cusparseHandle_t handle;
    cusparseCreate(&handle);

    int   *d_A_row_ptr, *d_A_col_ind;
    double *d_A_val;
    CUDA_CHECK(cudaMalloc(&d_A_row_ptr, sizeof(int) * (A.rows + 1)));
    CUDA_CHECK(cudaMalloc(&d_A_col_ind, sizeof(int) * A.nnz));
    CUDA_CHECK(cudaMalloc(&d_A_val,     sizeof(double) * A.nnz));
    CUDA_CHECK(cudaMemcpy(d_A_row_ptr, A.csr_row_ptr.data(),
                          sizeof(int) * (A.rows + 1), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_A_col_ind, A.csr_col_ind.data(),
                          sizeof(int) * A.nnz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_A_val, A.csr_val.data(),
                          sizeof(double) * A.nnz, cudaMemcpyHostToDevice));

    cusparseSpMatDescr_t desc_A;
    cusparseCreateCsr(&desc_A, A.rows, A.cols, A.nnz,
        d_A_row_ptr, d_A_col_ind, d_A_val,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);

    int *d_C_row_ptr;
    CUDA_CHECK(cudaMalloc(&d_C_row_ptr, sizeof(int) * (A.rows + 1)));

    // Algorithm list: user hint first, then all others
    cusparseSpGEMMAlg_t algs[] = {
        CUSPARSE_SPGEMM_DEFAULT,
        CUSPARSE_SPGEMM_ALG1,
        CUSPARSE_SPGEMM_ALG2,
        CUSPARSE_SPGEMM_ALG3
    };
    int n_algs = 4;
    // If user specified an algorithm, try it first
    if (alg_hint >= 0 && alg_hint <= 5) {
        cusparseSpGEMMAlg_t hint = (cusparseSpGEMMAlg_t)alg_hint;
        // Move hint to front
        for (int i = 0; i < n_algs; i++) {
            if (algs[i] == hint) { algs[i] = algs[0]; algs[0] = hint; break; }
        }
    }

    std::vector<double> times;
    int   *d_C_col_ind = nullptr;
    double *d_C_val = nullptr;
    int64_t last_C_nnz = 0;
    int alg_used = -1;

    for (int i = 0; i < n_algs; i++) {
        fprintf(stderr, "Trying alg=%d...\n", (int)algs[i]);
        if (d_C_col_ind) { cudaFree(d_C_col_ind); d_C_col_ind = nullptr; }
        if (d_C_val)     { cudaFree(d_C_val);     d_C_val = nullptr;     }

        bool ok = try_spgemm(handle, desc_A, A.rows, A.cols,
                             d_A_row_ptr, d_C_row_ptr,
                             1.0, 0.0, algs[i],
                             warmup, rounds,
                             (int&)last_C_nnz, d_C_col_ind, d_C_val, times);
        if (ok) {
            alg_used = (int)algs[i];
            fprintf(stderr, "  alg=%d succeeded! %d rounds, median %.2f ms\n",
                    alg_used, (int)times.size(),
                    times.empty() ? -1 : times[times.size()/2]);
            break;
        }
    }

    if (alg_used < 0 || times.empty()) {
        fprintf(stderr, "All algorithms failed!\n");
        res.success = false;
    } else {
        res.success = true;
        res.alg_used = alg_used;
        res.result_nnz = (int)last_C_nnz;
        res.nnz_per_row = (int)last_C_nnz / A.rows;

        std::sort(times.begin(), times.end());
        res.median_ms = times[times.size() / 2];
        res.min_ms = times.front();
        res.max_ms = times.back();

        // Download C for trace/threshold
        std::vector<int> h_C_row_ptr(A.rows + 1);
        CUDA_CHECK(cudaMemcpy(h_C_row_ptr.data(), d_C_row_ptr,
                              sizeof(int) * (A.rows + 1), cudaMemcpyDeviceToHost));
        std::vector<int> h_C_col_ind(last_C_nnz);
        std::vector<double> h_C_val(last_C_nnz);
        CUDA_CHECK(cudaMemcpy(h_C_col_ind.data(), d_C_col_ind,
                              sizeof(int) * last_C_nnz, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_C_val.data(), d_C_val,
                              sizeof(double) * last_C_nnz, cudaMemcpyDeviceToHost));

        double trace = 0.0, grand_sum = 0.0;
        int nnz_above = 0;
        for (int r = 0; r < A.rows; r++) {
            for (int j = h_C_row_ptr[r]; j < h_C_row_ptr[r + 1]; j++) {
                if (h_C_col_ind[j] == r) trace += h_C_val[j];
                grand_sum += h_C_val[j];
                if (std::abs(h_C_val[j]) > threshold) nnz_above++;
            }
        }
        res.trace = trace;
        res.grand_sum = grand_sum;
        res.nnz_above_threshold = nnz_above;

        fprintf(stderr, "Result nnz: %lld, nnz above threshold: %d, nnz/row: %d\n",
                (long long)last_C_nnz, nnz_above, res.nnz_per_row);
        fprintf(stderr, "Trace: %.15e, Grand sum: %.15e\n", trace, grand_sum);
    }

    // Cleanup
    if (d_C_col_ind) cudaFree(d_C_col_ind);
    if (d_C_val)     cudaFree(d_C_val);
    cusparseDestroySpMat(desc_A);
    cudaFree(d_A_row_ptr);
    cudaFree(d_A_col_ind);
    cudaFree(d_A_val);
    cudaFree(d_C_row_ptr);
    cusparseDestroy(handle);

    return res;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s matrix.bin [--rounds N] [--warmup N] [--threshold T] [--alg N]\n", argv[0]);
        return 1;
    }
    const char *matrix_file = argv[1];
    int   rounds    = 10;
    int   warmup    = 3;
    double threshold = 1e-4;
    int   alg_hint   = -1;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--rounds") == 0 && i + 1 < argc) rounds = atoi(argv[++i]);
        else if (strcmp(argv[i], "--warmup") == 0 && i + 1 < argc) warmup = atoi(argv[++i]);
        else if (strcmp(argv[i], "--threshold") == 0 && i + 1 < argc) threshold = atof(argv[++i]);
        else if (strcmp(argv[i], "--alg") == 0 && i + 1 < argc) alg_hint = atoi(argv[++i]);
    }

    int dev; cudaGetDevice(&dev);
    cudaDeviceProp prop; cudaGetDeviceProperties(&prop, dev);
    fprintf(stderr, "GPU: %s, SM %d.%d, %.0f MB, %d SMs\n",
            prop.name, prop.major, prop.minor,
            prop.totalGlobalMem / 1e6, prop.multiProcessorCount);

    COOMatrix coo = read_ntpoly_binary(matrix_file);
    CSRMatrix  A  = coo_to_csr(coo);
    fprintf(stderr, "CSR: %d x %d, nnz = %d\n", A.rows, A.cols, A.nnz);

    BenchResult res = benchmark_spgemm(A, warmup, rounds, threshold, alg_hint);

    printf("{\n");
    printf("  \"app\": \"bench_cusparse\",\n");
    printf("  \"gpu\": \"%s\",\n", prop.name);
    printf("  \"matrix\": \"%s\",\n", matrix_file);
    printf("  \"matrix_dimension\": %d,\n", A.rows);
    printf("  \"input_nnz\": %d,\n", A.nnz);
    printf("  \"input_nnz_per_row\": %.1f,\n", (double)A.nnz / A.rows);
    printf("  \"warmup\": %d,\n", warmup);
    printf("  \"rounds\": %d,\n", rounds);
    printf("  \"threshold\": %.15e,\n", threshold);
    if (res.success) {
        printf("  \"alg_used\": %d,\n", res.alg_used);
        printf("  \"gemm_median_ms\": %.4f,\n", res.median_ms);
        printf("  \"gemm_min_ms\": %.4f,\n", res.min_ms);
        printf("  \"gemm_max_ms\": %.4f,\n", res.max_ms);
        printf("  \"result_nnz\": %d,\n", res.result_nnz);
        printf("  \"result_nnz_per_row\": %d,\n", res.nnz_per_row);
        printf("  \"nnz_above_threshold\": %d,\n", res.nnz_above_threshold);
        printf("  \"trace\": %.15e,\n", res.trace);
        printf("  \"grand_sum\": %.15e\n", res.grand_sum);
    } else {
        printf("  \"success\": false,\n");
        printf("  \"error\": \"all algorithms failed\"\n");
    }
    printf("}\n");
    return 0;
}
