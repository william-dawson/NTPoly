// cusparse_gemm.cu  — optimized v4
// C wrapper for cuSPARSE SpGEMM with GPU-side thresholding.
//
// All GPU buffers are persistent across calls (no per-round cudaMalloc/cudaFree):
//   - A (uploaded once, reused while A doesn't change)
//   - C full-product (row_ptr, col_ind, val) — realloc only if C_nnz grows
//   - cuSPARSE work buffers (b1, b2) — realloc only if size grows
//   - Threshold pipeline buffers (keep, row_count, out_row_ptr, out_col_ind, out_val)
//
// Algorithm selection via env var CUSPARSE_ALG:
//   unset = try DEFAULT, then ALG1/ALG2/ALG3 (fallback chain)
//   0     = force DEFAULT
//   1     = force ALG1
//   2     = force ALG2
//   3     = force ALG3
//   -1    = benchmark ALL algorithms, print per-alg timing, use fastest

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <deque>
#include <cmath>
#include <ctime>
#include <utility>
#include <climits>
#include <cstdint>
#include <cusparse.h>
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>

#define CU_CHECK(call) do {                                              \
    cudaError_t err = (call);                                             \
    if (err != cudaSuccess) {                                             \
        fprintf(stderr, "CUDA error at %s:%d: %s\n",                      \
                __FILE__, __LINE__, cudaGetErrorString(err));              \
        return -1;                                                        \
    }                                                                     \
} while(0)

#define CS_CHECK(call) do {                                              \
    cusparseStatus_t s = (call);                                          \
    if (s != CUSPARSE_STATUS_SUCCESS) {                                   \
        fprintf(stderr, "cuSPARSE error at %s:%d: %s\n",                  \
                __FILE__, __LINE__, cusparseGetErrorString(s));            \
        return -1;                                                        \
    }                                                                     \
} while(0)

#define CS_BOOL(call) do { \
    cusparseStatus_t status = (call); \
    if (status != CUSPARSE_STATUS_SUCCESS) { \
        fprintf(stderr, "cuSPARSE error at %s:%d: %s\n", __FILE__, __LINE__, cusparseGetErrorString(status)); \
        return false; \
    } \
} while (0)

#ifdef CUSPARSE_PROFILE
static inline double now_sec() {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}
static double g_t0, g_t1;
// Rank tag for [prof] lines (mpirun sets OMPI_COMM_WORLD_RANK).
static const char *prof_rank() {
    static char buf[16];
    static bool init = false;
    if (!init) {
        const char *e = getenv("OMPI_COMM_WORLD_RANK");
        if (!e) e = getenv("PMI_RANK");
        if (!e) e = getenv("SLURM_PROCID");
        snprintf(buf, sizeof(buf), "%s", e ? e : "?");
        init = true;
    }
    return buf;
}
#define PROF_START() g_t0 = now_sec()
#define PROF_END(label) do { \
    g_t1 = now_sec(); \
    fprintf(stderr, "[prof r%s] %-30s %8.3f ms\n", prof_rank(), label, (g_t1 - g_t0) * 1e3); \
    g_t0 = g_t1; \
} while(0)
static cudaEvent_t g_ev0, g_ev1;
static bool g_ev_init = false;
static inline void ev_init() {
    if (!g_ev_init) { cudaEventCreate(&g_ev0); cudaEventCreate(&g_ev1); g_ev_init = true; }
}
#define GPU_START() do { ev_init(); cudaEventRecord(g_ev0); } while(0)
#define GPU_END(label) do { \
    cudaEventRecord(g_ev1); cudaEventSynchronize(g_ev1); \
    float ms = 0; cudaEventElapsedTime(&ms, g_ev0, g_ev1); \
    fprintf(stderr, "[prof r%s] %-30s %8.3f ms (GPU)\n", prof_rank(), label, ms); \
} while(0)
#define GPU_END_SILENT(ms_out) do { \
    cudaEventRecord(g_ev1); cudaEventSynchronize(g_ev1); \
    float _ms = 0; cudaEventElapsedTime(&_ms, g_ev0, g_ev1); \
    ms_out = _ms; \
} while(0)
#else
#define PROF_START()
#define PROF_END(label)
#define GPU_START()
#define GPU_END(label)
#define GPU_END_SILENT(ms_out)
#endif

// --- Persistent state (across calls) ---
static cusparseHandle_t g_handle = nullptr;
static bool g_handle_init = false;

// Persistent A on GPU
static int *d_A_row_ptr = nullptr, *d_A_col_ind = nullptr;
static double *d_A_val = nullptr;
static cusparseSpMatDescr_t g_desc_A = nullptr;
static int g_A_n = 0, g_A_nnz = 0;

// Persistent C (full product) buffers on GPU
static int *d_C_row_ptr = nullptr;
static int *d_C_col_ind = nullptr;
static double *d_C_val = nullptr;
static int64_t g_C_capacity = 0;   // capacity of d_C_col_ind / d_C_val
static int g_C_n = 0;              // matrix dimension for d_C_row_ptr
static cusparseSpMatDescr_t g_desc_C = nullptr;

// Persistent work buffers for cuSPARSE SpGEMM
static void *d_work1 = nullptr;
static size_t g_work1_sz = 0;
static void *d_work2 = nullptr;
static size_t g_work2_sz = 0;

// Warm flag: after first successful SpGEMM, skip size-query calls (A is same)
static bool g_spgemm_warm = false;

// --- Multi-GPU path: A_local * B (1D row decomposition) ---
// Persistent A_local (row-block of A) on GPU
static int *d_Aloc_row_ptr = nullptr, *d_Aloc_col_ind = nullptr;
static double *d_Aloc_val = nullptr;
static cusparseSpMatDescr_t g_desc_Aloc = nullptr;
static int g_Aloc_m = 0, g_Aloc_n = 0, g_Aloc_nnz = 0;

// Persistent B (full matrix) on GPU
static int *d_B_row_ptr = nullptr, *d_B_col_ind = nullptr;
static double *d_B_val = nullptr;
static cusparseSpMatDescr_t g_desc_B = nullptr;
static int g_B_n = 0, g_B_k = 0, g_B_nnz = 0;  // B is n × k

// Separate warm flag for A*B path
static bool g_spgemm_ab_warm = false;

// Persistent threshold-pipeline buffers on GPU
static int *d_keep = nullptr;
static int64_t g_keep_capacity = 0;
static int *d_row_count = nullptr;
static int g_row_count_n = 0;
static int *d_out_row_ptr = nullptr;
static int g_out_row_ptr_n = 0;

// Persistent output buffers on GPU
static int *d_out_col_ind = nullptr;
static double *d_out_val = nullptr;
static int g_out_capacity = 0;

// Persistent output buffers on host (avoid per-round malloc)
static int *h_out_row_ptr = nullptr;
static int *h_out_col_ind = nullptr;
static double *h_out_val = nullptr;
static int h_out_capacity_nnz = 0;
static int h_out_n = 0;

static bool init_gpu() {
    cudaGetLastError();
    // Select GPU based on local MPI rank (for multi-GPU scaling)
    int dev = 0;
    const char *lr = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    if (!lr) lr = getenv("SLURM_LOCALID");
    if (lr) dev = atoi(lr);
    int ndev = 0;
    cudaGetDeviceCount(&ndev);
    if (ndev > 0 && dev >= ndev) dev = dev % ndev;
    if (cudaSetDevice(dev) != cudaSuccess) return false;
    int cur = -1;
    if (cudaGetDevice(&cur) != cudaSuccess || cur != dev) return false;
    return true;
}

static cusparseHandle_t get_handle() {
    if (!g_handle_init) {
        if (!init_gpu()) { g_handle = nullptr; g_handle_init = true; return nullptr; }
        if (cusparseCreate(&g_handle) != CUSPARSE_STATUS_SUCCESS) g_handle = nullptr;
        g_handle_init = true;
    }
    return g_handle;
}

// --- GPU kernels ---

// Subtract 1 from each int element (1-based -> 0-based conversion on GPU).
__global__ void subtract_one_kernel(int *data, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) data[idx] -= 1;
}

__global__ void threshold_mark_kernel(
    const double *val, int nnz, double threshold, int *keep)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nnz)
        keep[idx] = (fabs(val[idx]) > threshold) ? 1 : 0;
}

__global__ void count_per_row_kernel(
    const int *row_ptr, const int *keep, int n, int *row_count)
{
    int row = blockIdx.x;
    if (row >= n) return;
    int start = row_ptr[row];
    int end = row_ptr[row + 1];
    int tid = threadIdx.x;
    int local = 0;
    for (int j = start + tid; j < end; j += blockDim.x)
        local += keep[j];
    __shared__ int sdata[256];
    sdata[tid] = local;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    if (tid == 0) row_count[row] = sdata[0];
}

__global__ void compact_kernel(
    const int *row_ptr_in, const int *col_ind_in, const double *val_in,
    const int *keep, const int *row_ptr_out, int n,
    int *col_ind_out, double *val_out)
{
    int row = blockIdx.x;
    if (row >= n) return;
    int out_pos = row_ptr_out[row];
    int start = row_ptr_in[row];
    int end = row_ptr_in[row + 1];
    for (int j = start; j < end; j++) {
        if (keep[j]) {
            col_ind_out[out_pos] = col_ind_in[j] + 1;  // 0-based -> 1-based (NTPoly)
            val_out[out_pos] = val_in[j];
            out_pos++;
        }
    }
}

// Upload A to GPU if not already uploaded (or if size changed).
static int ensure_a_uploaded(int n, int nnz_a,
    const int *row_ptr_a, const int *col_ind_a, const double *val_a)
{
    if (d_A_row_ptr && g_A_n == n && g_A_nnz == nnz_a)
        return 0;  // already uploaded

    // A changed — reset warm flag (work buffer sizes may differ)
    g_spgemm_warm = false;

    if (d_A_row_ptr) { cudaFree(d_A_row_ptr); d_A_row_ptr = nullptr; }
    if (d_A_col_ind) { cudaFree(d_A_col_ind); d_A_col_ind = nullptr; }
    if (d_A_val) { cudaFree(d_A_val); d_A_val = nullptr; }
    if (g_desc_A) { cusparseDestroySpMat(g_desc_A); g_desc_A = nullptr; }

    CU_CHECK(cudaMalloc(&d_A_row_ptr, sizeof(int) * (n + 1)));
    CU_CHECK(cudaMalloc(&d_A_col_ind, sizeof(int) * nnz_a));
    CU_CHECK(cudaMalloc(&d_A_val, sizeof(double) * nnz_a));

    CU_CHECK(cudaMemcpy(d_A_row_ptr, row_ptr_a, sizeof(int) * (n + 1),
                        cudaMemcpyHostToDevice));
    CU_CHECK(cudaMemcpy(d_A_val, val_a, sizeof(double) * nnz_a,
                        cudaMemcpyHostToDevice));
    // NTPoly inner_index is 1-based; cuSPARSE INDEX_BASE_ZERO needs 0-based.
    // Copy then subtract 1 on GPU (runs once, not on critical path).
    CU_CHECK(cudaMemcpy(d_A_col_ind, col_ind_a, sizeof(int) * nnz_a,
                        cudaMemcpyHostToDevice));
    {
        int t = 256, b = (nnz_a + t - 1) / t;
        subtract_one_kernel<<<b, t>>>(d_A_col_ind, nnz_a);
        CU_CHECK(cudaGetLastError());
        CU_CHECK(cudaDeviceSynchronize());
    }

    CS_CHECK(cusparseCreateCsr(&g_desc_A, n, n, nnz_a,
        d_A_row_ptr, d_A_col_ind, d_A_val,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

    g_A_n = n;
    g_A_nnz = nnz_a;
    return 0;
}

// Reuse A_local allocations, but always upload the supplied matrix.
// A_local is m × n (m rows, n columns), CSR format.
static int ensure_aloc_uploaded(int m, int n, int nnz_a,
    const int *row_ptr_a, const int *col_ind_a, const double *val_a)
{
    if (d_Aloc_row_ptr && g_Aloc_m == m && g_Aloc_n == n && g_Aloc_nnz == nnz_a) {
        g_spgemm_ab_warm = false;
        CU_CHECK(cudaMemcpy(d_Aloc_row_ptr, row_ptr_a, sizeof(int)*(size_t)(m+1), cudaMemcpyHostToDevice));
        if (nnz_a > 0) {
            CU_CHECK(cudaMemcpy(d_Aloc_col_ind, col_ind_a, sizeof(int)*(size_t)nnz_a, cudaMemcpyHostToDevice));
            CU_CHECK(cudaMemcpy(d_Aloc_val, val_a, sizeof(double)*(size_t)nnz_a, cudaMemcpyHostToDevice));
            subtract_one_kernel<<<(nnz_a+255)/256, 256>>>(d_Aloc_col_ind, nnz_a);
            CU_CHECK(cudaGetLastError());
            CU_CHECK(cudaDeviceSynchronize());
        }
        return 0;
    }

    g_spgemm_ab_warm = false;

    if (d_Aloc_row_ptr) { cudaFree(d_Aloc_row_ptr); d_Aloc_row_ptr = nullptr; }
    if (d_Aloc_col_ind) { cudaFree(d_Aloc_col_ind); d_Aloc_col_ind = nullptr; }
    if (d_Aloc_val) { cudaFree(d_Aloc_val); d_Aloc_val = nullptr; }
    if (g_desc_Aloc) { cusparseDestroySpMat(g_desc_Aloc); g_desc_Aloc = nullptr; }

    CU_CHECK(cudaMalloc(&d_Aloc_row_ptr, sizeof(int) * (m + 1)));
    CU_CHECK(cudaMalloc(&d_Aloc_col_ind, sizeof(int) * nnz_a));
    CU_CHECK(cudaMalloc(&d_Aloc_val, sizeof(double) * nnz_a));

    CU_CHECK(cudaMemcpy(d_Aloc_row_ptr, row_ptr_a, sizeof(int) * (m + 1),
                        cudaMemcpyHostToDevice));
    CU_CHECK(cudaMemcpy(d_Aloc_val, val_a, sizeof(double) * nnz_a,
                        cudaMemcpyHostToDevice));
    CU_CHECK(cudaMemcpy(d_Aloc_col_ind, col_ind_a, sizeof(int) * nnz_a,
                        cudaMemcpyHostToDevice));
    {
        int t = 256, b = (nnz_a + t - 1) / t;
        subtract_one_kernel<<<b, t>>>(d_Aloc_col_ind, nnz_a);
        CU_CHECK(cudaGetLastError());
        CU_CHECK(cudaDeviceSynchronize());
    }

    CS_CHECK(cusparseCreateCsr(&g_desc_Aloc, m, n, nnz_a,
        d_Aloc_row_ptr, d_Aloc_col_ind, d_Aloc_val,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

    g_Aloc_n = n;
    g_Aloc_m = m;
    g_Aloc_nnz = nnz_a;
    return 0;
}

// Reuse B allocations, but always upload the supplied matrix.
// B is n × k, CSR format.
static int ensure_b_uploaded(int n, int k, int nnz_b,
    const int *row_ptr_b, const int *col_ind_b, const double *val_b)
{
    if (d_B_row_ptr && g_B_n == n && g_B_k == k && g_B_nnz == nnz_b) {
        g_spgemm_ab_warm = false;
        CU_CHECK(cudaMemcpy(d_B_row_ptr, row_ptr_b, sizeof(int)*(size_t)(n+1), cudaMemcpyHostToDevice));
        if (nnz_b > 0) {
            CU_CHECK(cudaMemcpy(d_B_col_ind, col_ind_b, sizeof(int)*(size_t)nnz_b, cudaMemcpyHostToDevice));
            CU_CHECK(cudaMemcpy(d_B_val, val_b, sizeof(double)*(size_t)nnz_b, cudaMemcpyHostToDevice));
            subtract_one_kernel<<<(nnz_b+255)/256, 256>>>(d_B_col_ind, nnz_b);
            CU_CHECK(cudaGetLastError());
            CU_CHECK(cudaDeviceSynchronize());
        }
        return 0;
    }

    g_spgemm_ab_warm = false;

    if (d_B_row_ptr) { cudaFree(d_B_row_ptr); d_B_row_ptr = nullptr; }
    if (d_B_col_ind) { cudaFree(d_B_col_ind); d_B_col_ind = nullptr; }
    if (d_B_val) { cudaFree(d_B_val); d_B_val = nullptr; }
    if (g_desc_B) { cusparseDestroySpMat(g_desc_B); g_desc_B = nullptr; }

    CU_CHECK(cudaMalloc(&d_B_row_ptr, sizeof(int) * (n + 1)));
    CU_CHECK(cudaMalloc(&d_B_col_ind, sizeof(int) * nnz_b));
    CU_CHECK(cudaMalloc(&d_B_val, sizeof(double) * nnz_b));

    CU_CHECK(cudaMemcpy(d_B_row_ptr, row_ptr_b, sizeof(int) * (n + 1),
                        cudaMemcpyHostToDevice));
    CU_CHECK(cudaMemcpy(d_B_val, val_b, sizeof(double) * nnz_b,
                        cudaMemcpyHostToDevice));
    CU_CHECK(cudaMemcpy(d_B_col_ind, col_ind_b, sizeof(int) * nnz_b,
                        cudaMemcpyHostToDevice));
    {
        int t = 256, b = (nnz_b + t - 1) / t;
        subtract_one_kernel<<<b, t>>>(d_B_col_ind, nnz_b);
        CU_CHECK(cudaGetLastError());
        CU_CHECK(cudaDeviceSynchronize());
    }

    CS_CHECK(cusparseCreateCsr(&g_desc_B, n, k, nnz_b,
        d_B_row_ptr, d_B_col_ind, d_B_val,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

    g_B_n = n;
    g_B_k = k;
    g_B_nnz = nnz_b;
    return 0;
}

// Ensure persistent C row_ptr is allocated for dimension n.
static int ensure_c_row_ptr(int n) {
    if (d_C_row_ptr && g_C_n == n) return 0;
    if (d_C_row_ptr) { cudaFree(d_C_row_ptr); d_C_row_ptr = nullptr; }
    if (g_desc_C) { cusparseDestroySpMat(g_desc_C); g_desc_C = nullptr; }
    CU_CHECK(cudaMalloc(&d_C_row_ptr, sizeof(int) * (n + 1)));
    g_C_n = n;
    return 0;
}

// Ensure persistent C col_ind/val have capacity >= needed.
static int ensure_c_data(int64_t needed) {
    if (g_C_capacity >= needed) return 0;
    if (d_C_col_ind) { cudaFree(d_C_col_ind); d_C_col_ind = nullptr; }
    if (d_C_val) { cudaFree(d_C_val); d_C_val = nullptr; }
    CU_CHECK(cudaMalloc(&d_C_col_ind, sizeof(int) * needed));
    CU_CHECK(cudaMalloc(&d_C_val, sizeof(double) * needed));
    g_C_capacity = needed;
    return 0;
}

// Ensure persistent work buffer has capacity >= needed.
static int ensure_work_buf(void **buf, size_t *cur_sz, size_t needed) {
    if (*cur_sz >= needed) return 0;
    if (*buf) { cudaFree(*buf); *buf = nullptr; }
    CU_CHECK(cudaMalloc(buf, needed));
    *cur_sz = needed;
    return 0;
}

// Ensure persistent keep buffer has capacity >= needed.
static int ensure_keep(int64_t needed) {
    if (g_keep_capacity >= needed) return 0;
    if (d_keep) { cudaFree(d_keep); d_keep = nullptr; }
    CU_CHECK(cudaMalloc(&d_keep, sizeof(int) * needed));
    g_keep_capacity = needed;
    return 0;
}

// Ensure persistent row_count buffer for dimension n.
static int ensure_row_count(int n) {
    if (d_row_count && g_row_count_n == n) return 0;
    if (d_row_count) { cudaFree(d_row_count); d_row_count = nullptr; }
    CU_CHECK(cudaMalloc(&d_row_count, sizeof(int) * n));
    g_row_count_n = n;
    return 0;
}

// Ensure persistent out_row_ptr for dimension n.
static int ensure_out_row_ptr(int n) {
    if (d_out_row_ptr && g_out_row_ptr_n == n) return 0;
    if (d_out_row_ptr) { cudaFree(d_out_row_ptr); d_out_row_ptr = nullptr; }
    CU_CHECK(cudaMalloc(&d_out_row_ptr, sizeof(int) * (n + 1)));
    g_out_row_ptr_n = n;
    return 0;
}

// Ensure persistent out col_ind/val have capacity >= needed.
static int ensure_out_data(int needed) {
    if (g_out_capacity >= needed) return 0;
    if (d_out_col_ind) { cudaFree(d_out_col_ind); d_out_col_ind = nullptr; }
    if (d_out_val) { cudaFree(d_out_val); d_out_val = nullptr; }
    CU_CHECK(cudaMalloc(&d_out_col_ind, sizeof(int) * needed));
    CU_CHECK(cudaMalloc(&d_out_val, sizeof(double) * needed));
    g_out_capacity = needed;
    return 0;
}

// Run one SpGEMM algorithm with two arbitrary descriptors.
// C = desc_A * desc_B, where C is c_rows × c_cols.
// Returns true on success, fills C_nnz.
static bool run_spgemm_general(cusparseHandle_t handle,
    int c_rows, int c_cols,
    cusparseSpMatDescr_t desc_A, cusparseSpMatDescr_t desc_B,
    cusparseSpGEMMAlg_t alg, int64_t *C_nnz_out,
    bool *warm_flag, bool accept_empty = false)
{
    double alpha = 1.0, beta = 0.0;

#ifdef CUSPARSE_PROFILE
    static cudaEvent_t pe0, pe1;
    static bool pe_init = false;
    if (!pe_init) {
        cudaEventCreate(&pe0); cudaEventCreate(&pe1);
        pe_init = true;
    }
    #define PHASE_GPU_START() cudaEventRecord(pe0)
    #define PHASE_GPU_END(label) do { \
        cudaEventRecord(pe1); cudaEventSynchronize(pe1); \
        float _ms = 0; cudaEventElapsedTime(&_ms, pe0, pe1); \
        fprintf(stderr, "[prof r%s]   %-28s %8.3f ms (GPU)\n", prof_rank(), label, _ms); \
    } while(0)
#else
    #define PHASE_GPU_START()
    #define PHASE_GPU_END(label)
#endif

    // Recreate desc_C (nnz starts at 0 each call)
    if (g_desc_C) cusparseDestroySpMat(g_desc_C);
    CS_BOOL(cusparseCreateCsr(&g_desc_C, c_rows, c_cols, 0,
        d_C_row_ptr, nullptr, nullptr,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

    cusparseSpGEMMDescr_t desc;
    CS_BOOL(cusparseSpGEMM_createDescr(&desc));

    // Phase 1: work estimation
    if (*warm_flag && g_work1_sz > 0) {
        size_t b1sz = g_work1_sz;
        if (cusparseSpGEMM_workEstimation(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc_A, desc_B, &beta,
            g_desc_C, CUDA_R_64F, alg, desc, &b1sz, d_work1)
            != CUSPARSE_STATUS_SUCCESS) {
            cusparseSpGEMM_destroyDescr(desc);
            *warm_flag = false;
            return false;
        }
    } else {
        size_t b1sz = 0;
        if (cusparseSpGEMM_workEstimation(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc_A, desc_B, &beta,
            g_desc_C, CUDA_R_64F, alg, desc, &b1sz, nullptr)
            != CUSPARSE_STATUS_SUCCESS) {
            cusparseSpGEMM_destroyDescr(desc);
            return false;
        }
        if (ensure_work_buf(&d_work1, &g_work1_sz, b1sz) != 0) {
            cusparseSpGEMM_destroyDescr(desc);
            return false;
        }
        if (cusparseSpGEMM_workEstimation(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc_A, desc_B, &beta,
            g_desc_C, CUDA_R_64F, alg, desc, &b1sz, d_work1)
            != CUSPARSE_STATUS_SUCCESS) {
            cusparseSpGEMM_destroyDescr(desc);
            return false;
        }
    }

    // Phase 2: compute
    if (*warm_flag && g_work2_sz > 0) {
        size_t b2sz = g_work2_sz;
        PHASE_GPU_START();
        if (cusparseSpGEMM_compute(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc_A, desc_B, &beta,
            g_desc_C, CUDA_R_64F, alg, desc, &b2sz, d_work2)
            != CUSPARSE_STATUS_SUCCESS) {
            cusparseSpGEMM_destroyDescr(desc);
            *warm_flag = false;
            return false;
        }
        if (cudaDeviceSynchronize() != cudaSuccess) return false;
        PHASE_GPU_END("  compute (warm)");
    } else {
        size_t b2sz = 0;
        if (cusparseSpGEMM_compute(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc_A, desc_B, &beta,
            g_desc_C, CUDA_R_64F, alg, desc, &b2sz, nullptr)
            != CUSPARSE_STATUS_SUCCESS) {
            cusparseSpGEMM_destroyDescr(desc);
            return false;
        }
        if (ensure_work_buf(&d_work2, &g_work2_sz, b2sz) != 0) {
            cusparseSpGEMM_destroyDescr(desc);
            return false;
        }
        PHASE_GPU_START();
        if (cusparseSpGEMM_compute(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc_A, desc_B, &beta,
            g_desc_C, CUDA_R_64F, alg, desc, &b2sz, d_work2)
            != CUSPARSE_STATUS_SUCCESS) {
            cusparseSpGEMM_destroyDescr(desc);
            return false;
        }
        if (cudaDeviceSynchronize() != cudaSuccess) return false;
        PHASE_GPU_END("  compute (cold)");
    }

    // Get C_nnz
    int64_t Cr, Cc, C_nnz;
    CS_BOOL(cusparseSpMatGetSize(g_desc_C, &Cr, &Cc, &C_nnz));
    if (C_nnz == 0) {
        cusparseSpGEMM_destroyDescr(desc);
        *C_nnz_out = 0;
        return accept_empty;
    }
    if (C_nnz > INT_MAX) {
        cusparseSpGEMM_destroyDescr(desc);
        return false;
    }

    // Ensure C data buffers are large enough
    if (ensure_c_data(C_nnz) != 0) {
        cusparseSpGEMM_destroyDescr(desc);
        return false;
    }
    CS_BOOL(cusparseCsrSetPointers(g_desc_C, d_C_row_ptr, d_C_col_ind, d_C_val));

    // Phase 3: copy
    PHASE_GPU_START();
    if (cusparseSpGEMM_copy(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc_A, desc_B, &beta,
        g_desc_C, CUDA_R_64F, alg, desc)
        != CUSPARSE_STATUS_SUCCESS) {
        cusparseSpGEMM_destroyDescr(desc);
        return false;
    }
    if (cudaDeviceSynchronize() != cudaSuccess) return false;
    PHASE_GPU_END("  copy");

    cusparseSpGEMM_destroyDescr(desc);
    *C_nnz_out = C_nnz;
    *warm_flag = true;
    return true;
}

// Wrapper for A*A path (calls run_spgemm_general with desc_A for both operands)
static bool run_spgemm_one(cusparseHandle_t handle, int n,
    cusparseSpGEMMAlg_t alg, int64_t *C_nnz_out)
{
    return run_spgemm_general(handle, n, n, g_desc_A, g_desc_A,
        alg, C_nnz_out, &g_spgemm_warm);
}

// Run SpGEMM with algorithm selection (env var CUSPARSE_ALG).
// Tries DEFAULT→ALG1→ALG2→ALG3 by default, or forces one, or benchmarks all.
static bool do_spgemm(cusparseHandle_t handle,
    int c_rows, int c_cols,
    cusparseSpMatDescr_t desc_A, cusparseSpMatDescr_t desc_B,
    bool *warm_flag,
    int64_t *C_nnz_out)
{
    static int alg_select = -999;
    if (alg_select == -999) {
        const char *env = getenv("CUSPARSE_ALG");
        alg_select = env ? atoi(env) : -1;
    }

    cusparseSpGEMMAlg_t algs[] = {
        CUSPARSE_SPGEMM_DEFAULT, CUSPARSE_SPGEMM_ALG1,
        CUSPARSE_SPGEMM_ALG2, CUSPARSE_SPGEMM_ALG3
    };
    const char *alg_names[] = {"DEFAULT", "ALG1", "ALG2", "ALG3"};

    bool success = false;

    GPU_START();

    if (alg_select == -1) {
        for (int i = 0; i < 4 && !success; i++)
            success = run_spgemm_general(handle, c_rows, c_cols,
                desc_A, desc_B, algs[i], C_nnz_out, warm_flag);
    } else if (alg_select >= 0 && alg_select <= 3) {
        success = run_spgemm_general(handle, c_rows, c_cols,
            desc_A, desc_B, algs[alg_select], C_nnz_out, warm_flag);
    } else if (alg_select == -2) {
        int64_t best_nnz = 0;
        float best_ms = 1e30;
        int best_alg = -1;
        for (int i = 0; i < 4; i++) {
            cudaEvent_t ev_a, ev_b;
            cudaEventCreate(&ev_a); cudaEventCreate(&ev_b);
            cudaEventRecord(ev_a);
            int64_t this_nnz = 0;
            bool ok = run_spgemm_general(handle, c_rows, c_cols,
                desc_A, desc_B, algs[i], &this_nnz, warm_flag);
            cudaEventRecord(ev_b);
            cudaEventSynchronize(ev_b);
            float ms = 0;
            cudaEventElapsedTime(&ms, ev_a, ev_b);
            cudaEventDestroy(ev_a); cudaEventDestroy(ev_b);
            if (ok) {
                fprintf(stderr, "[alg-bench] %-10s  %8.3f ms  C_nnz=%lld  %s\n",
                    alg_names[i], ms, (long long)this_nnz,
                    ms < best_ms ? " <-- BEST" : "");
                if (ms < best_ms) { best_ms = ms; best_alg = i; best_nnz = this_nnz; }
            } else {
                fprintf(stderr, "[alg-bench] %-10s  FAILED\n", alg_names[i]);
            }
        }
        if (best_alg >= 0)
            success = run_spgemm_general(handle, c_rows, c_cols,
                desc_A, desc_B, algs[best_alg], C_nnz_out, warm_flag);
    }

    GPU_END("cuSPARSE SpGEMM");
    return success;
}

// Threshold pipeline: mark + count + scan + compact + download.
// Operates on d_C_row_ptr/d_C_col_ind/d_C_val (the full product on GPU).
// c_rows = number of rows of C (for kernel launch and scan).
static int do_threshold_pipeline(
    int c_rows, int64_t C_nnz, double threshold,
    int *nnz_c, int **row_ptr_c_out, int **col_ind_c_out, double **val_c_out)
{
    PROF_START();

    if (ensure_keep(C_nnz) != 0) return -1;
    if (ensure_row_count(c_rows) != 0) return -1;
    if (ensure_out_row_ptr(c_rows) != 0) return -1;

    {
        int t = 256, b = (C_nnz + t - 1) / t;
        threshold_mark_kernel<<<b, t>>>(d_C_val, (int)C_nnz, threshold, d_keep);
        CU_CHECK(cudaGetLastError());
    }

    count_per_row_kernel<<<c_rows, 256>>>(d_C_row_ptr, d_keep, c_rows, d_row_count);
    CU_CHECK(cudaGetLastError());
    CU_CHECK(cudaDeviceSynchronize());

    PROF_END("GPU mark + count");

    PROF_START();
    std::vector<int> h_row_count(c_rows);
    CU_CHECK(cudaMemcpy(h_row_count.data(), d_row_count,
                        sizeof(int) * c_rows, cudaMemcpyDeviceToHost));

    if (h_out_n < c_rows + 1) {
        if (h_out_row_ptr) free(h_out_row_ptr);
        h_out_row_ptr = (int *)malloc(sizeof(int) * (c_rows + 1));
        h_out_n = c_rows + 1;
    }
    h_out_row_ptr[0] = 0;
    int total_keep = 0;
    for (int r = 0; r < c_rows; r++) {
        total_keep += h_row_count[r];
        h_out_row_ptr[r + 1] = total_keep;
    }
    *nnz_c = total_keep;

    PROF_END("CPU scan + row_ptr");

    PROF_START();
    CU_CHECK(cudaMemcpy(d_out_row_ptr, h_out_row_ptr,
                        sizeof(int) * (c_rows + 1), cudaMemcpyHostToDevice));
    if (ensure_out_data(total_keep) != 0) return -1;

    compact_kernel<<<c_rows, 1>>>(
        d_C_row_ptr, d_C_col_ind, d_C_val, d_keep, d_out_row_ptr, c_rows,
        d_out_col_ind, d_out_val);
    CU_CHECK(cudaGetLastError());
    CU_CHECK(cudaDeviceSynchronize());

    PROF_END("GPU compact");

    PROF_START();
    if (h_out_capacity_nnz < total_keep) {
        if (h_out_col_ind) free(h_out_col_ind);
        if (h_out_val) free(h_out_val);
        h_out_col_ind = (int *)malloc(sizeof(int) * total_keep);
        h_out_val = (double *)malloc(sizeof(double) * total_keep);
        h_out_capacity_nnz = total_keep;
    }
    CU_CHECK(cudaMemcpy(h_out_col_ind, d_out_col_ind,
                        sizeof(int) * total_keep, cudaMemcpyDeviceToHost));
    CU_CHECK(cudaMemcpy(h_out_val, d_out_val,
                        sizeof(double) * total_keep, cudaMemcpyDeviceToHost));
    *row_ptr_c_out = h_out_row_ptr;
    *col_ind_c_out = h_out_col_ind;
    *val_c_out = h_out_val;

    PROF_END("D2H download (thresholded)");
    return 0;
}

// --- 2.5D k-split path: merge of ps partial products ---
// Partial products (one per slice, CSR) are concatenated and offset-fixed
// on the host, uploaded, scatter-added into a dense m x kk buffer, then
// thresholded and compacted back to CSR.  The dense buffer is persistent;
// it is sized m*kk*8 bytes, which only fits HBM for the small-matrix
// regime this path targets (the Fortran side guards the size too).

// Persistent k-split device buffers
static int *d_ks_rp = nullptr;         // ps*(m+1) row_ptr entries
static int g_ks_rp_cap = 0;
static int *d_ks_ci = nullptr;         // total_nnz col indices
static double *d_ks_v = nullptr;       // total_nnz values
static int64_t g_ks_cap = 0;
static double *d_ks_dense = nullptr;   // m*kk dense accumulator
static long long g_ks_dense_n = -1;
static unsigned *d_solver_bitmap=nullptr;
static size_t g_solver_bitmap_capacity=0;
static int g_solver_bitmap_rows=0,g_solver_bitmap_columns=0;
static bool g_solver_bitmap_valid=false;
// Solver-only streaming accumulation: consume each device partial before
// cuSPARSE reuses its output buffers. Legacy k-split callers retain host output.
static bool g_solver_accumulating = false;
static int g_solver_rows = 0, g_solver_cols = 0;

// Persistent host buffers for the partial product download
static int *h_p_row_ptr = nullptr;
static int h_p_row_n = 0;
static int *h_p_col_ind = nullptr;
static double *h_p_val = nullptr;
static int h_p_cap = 0;

// One block per (slice, row): grid (m, ps).  A (row,col) pair is unique
// within one slice's partial, but can appear in several slices -> atomicAdd.
// NOTE: col_ind_all is the RAW cuSPARSE output, 0-based (the 1-based NTPoly
// conversion normally done by compact_kernel does not apply to the
// unthresholded partial product).
__global__ void ksplit_scatter_kernel(
    int m, int kk,
    const int *row_ptr_all, const int *col_ind_all, const double *val_all,
    double *dense)
{
    int r = blockIdx.x;
    int s = blockIdx.y;
    const int *rp = row_ptr_all + (size_t)s * (m + 1);
    int start = rp[r];
    int end = rp[r + 1];
    for (int j = start + threadIdx.x; j < end; j += blockDim.x) {
        long long idx = (long long)r * kk + col_ind_all[j];
        atomicAdd(&dense[idx], val_all[j]);
    }
}

__global__ void ksplit_count_kernel(
    const double *dense, int m, int kk, double threshold, int *row_count)
{
    int row = blockIdx.x;
    if (row >= m) return;
    int tid = threadIdx.x;
    int local = 0;
    for (int c = tid; c < kk; c += blockDim.x)
        if (fabs(dense[(long long)row * kk + c]) > threshold) local++;
    __shared__ int sdata[256];
    sdata[tid] = local;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    if (tid == 0) row_count[row] = sdata[0];
}

// A single cuSPARSE partial has unique columns within each row. Chunks
// execute in order on the default stream, so no atomic operation is needed.
__global__ void solver_accumulate_kernel(int m, int k, const int *rp,
    const int *ci, const double *values, double *dense)
{
    int row = blockIdx.x;
    for (int j = rp[row] + threadIdx.x; j < rp[row + 1]; j += blockDim.x)
        dense[(long long)row * k + ci[j]] += values[j];
}

// Direct sparse multiply into a dense row accumulator. Four warps per row
// visit distinct A entries and coalesce reads along each corresponding B row.
// Different A entries may contribute to the same C entry, hence atomicAdd.
__global__ void solver_dense_product_kernel(int k, const int *arp,
    const int *aci, const double *av, const int *brp, const int *bci,
    const double *bv, double *dense)
{
    const int row = blockIdx.x, lane = threadIdx.x & 31;
    const int warp = threadIdx.x >> 5;
    for (int a = arp[row] + warp; a < arp[row + 1]; a += blockDim.x / 32) {
        const int shared = aci[a];
        const double value = av[a];
        for (int b = brp[shared] + lane; b < brp[shared + 1]; b += 32)
            atomicAdd(&dense[(long long)row * k + bci[b]], value * bv[b]);
    }
}

// A touched-entry bitmap lets subsequent passes skip the untouched dense
// accumulator. Payload and bitmap updates complete before count/compaction.
__global__ void solver_bitmap_product_kernel(int k,int words,const int *arp,
    const int *aci,const double *av,const int *brp,const int *bci,
    const double *bv,double *dense,unsigned *bitmap)
{
    int row=blockIdx.x,lane=threadIdx.x&31,warp=threadIdx.x>>5;
    for(int a=arp[row]+warp;a<arp[row+1];a+=blockDim.x/32) {
        int shared=aci[a]; double value=av[a];
        for(int b=brp[shared]+lane;b<brp[shared+1];b+=32) {
            int col=bci[b];
            atomicAdd(&dense[(long long)row*k+col],value*bv[b]);
            atomicOr(&bitmap[size_t(row)*words+col/32],1u<<(col&31));
        }
    }
}
__global__ void solver_bitmap_clear_kernel(int k,int words,double *dense,unsigned *bitmap)
{
    int row=blockIdx.x,lane=threadIdx.x&31;
    for(int word=threadIdx.x/32;word<words;word+=blockDim.x/32) {
        unsigned mask=lane==0?bitmap[size_t(row)*words+word]:0;
        mask=__shfl_sync(0xffffffffu,mask,0);
        int col=word*32+lane;
        if(col<k && (mask&(1u<<lane))) dense[(long long)row*k+col]=0;
        if(lane==0) bitmap[size_t(row)*words+word]=0;
    }
}
__global__ void solver_bitmap_count_kernel(int k,int words,const double *dense,
    const unsigned *bitmap,double threshold,int *counts)
{
    int row=blockIdx.x,tid=threadIdx.x,count=0;
    for(int word=tid;word<words;word+=blockDim.x) {
        unsigned mask=bitmap[size_t(row)*words+word];
        while(mask) {
            int col=word*32+__ffs(mask)-1;
            if(col<k && fabs(dense[(long long)row*k+col])>threshold) ++count;
            mask&=mask-1;
        }
    }
    __shared__ int scratch[256]; scratch[tid]=count; __syncthreads();
    for(int step=128;step;step/=2) {
        if(tid<step) scratch[tid]+=scratch[tid+step];
        __syncthreads();
    }
    if(tid==0) counts[row]=scratch[0];
}
__global__ void solver_bitmap_compact_kernel(int k,int words,const double *dense,
    const unsigned *bitmap,double threshold,const int *rp,int *ci,double *values)
{
    int row=blockIdx.x,lane=threadIdx.x,output=rp[row];
    for(int word=0;word<words;++word) {
        unsigned touched=lane==0?bitmap[size_t(row)*words+word]:0;
        touched=__shfl_sync(0xffffffffu,touched,0);
        if(!touched) continue;
        int col=word*32+lane;
        double value=col<k && (touched&(1u<<lane))?dense[(long long)row*k+col]:0;
        bool keep=col<k && (touched&(1u<<lane)) && fabs(value)>threshold;
        unsigned mask=__ballot_sync(0xffffffffu,keep);
        int offset=__popc(mask&((1u<<lane)-1));
        if(keep) { ci[output+offset]=col+1; values[output+offset]=value; }
        output+=__popc(mask);
    }
}

__global__ void ksplit_compact_kernel(
    const double *dense, int m, int kk, double threshold,
    const int *row_ptr_out, int *col_ind_out, double *val_out)
{
    int row = blockIdx.x;
    if (row >= m) return;
    // One complete warp per row: coalesced reads and ordered compaction.
    int lane = threadIdx.x;
    int out_pos = row_ptr_out[row];
    for (int tile = 0; tile < kk; tile += 32) {
        int c = tile + lane;
        double v = c < kk ? dense[(long long)row * kk + c] : 0.0;
        bool keep = c < kk && fabs(v) > threshold;
        unsigned mask = __ballot_sync(0xffffffffu, keep);
        int offset = __popc(mask & ((1u << lane) - 1));
        if (keep) {
            col_ind_out[out_pos + offset] = c + 1;
            val_out[out_pos + offset] = v;
        }
        out_pos += __popc(mask);
    }
}

extern "C" {

// Nsight can capture only the new solver loop, excluding the reference.
void gpu_profile_begin_() { (void)cudaProfilerStart(); }
void gpu_profile_end_() { (void)cudaProfilerStop(); }

int cusparse_solver_begin_(const int m, const int k) {
    if (g_solver_accumulating || m < 1 || k < 1) return -1;
    g_solver_bitmap_valid=false;
    const long long elements = (long long)m * k;
    if (elements > (long long)12e9) return -1;
    if (g_ks_dense_n != elements) {
        if (d_ks_dense) { CU_CHECK(cudaFree(d_ks_dense)); d_ks_dense = nullptr; g_ks_dense_n = -1; }
        CU_CHECK(cudaMalloc(&d_ks_dense, sizeof(double) * (size_t)elements));
        g_ks_dense_n = elements;
    }
    CU_CHECK(cudaMemset(d_ks_dense, 0, sizeof(double) * (size_t)elements));
    g_solver_rows = m; g_solver_cols = k;
    g_solver_accumulating = true;
    return 0;
}

int cusparse_gemm_threshold_(
    const int n, const int nnz_a,
    const int *row_ptr_a,   // 0-based (NTPoly outer_index)
    const int *col_ind_a,   // 1-based (NTPoly inner_index)
    const double *val_a,
    const double threshold,
    int *nnz_c,
    int **row_ptr_c_out, int **col_ind_c_out, double **val_c_out)
{
    cusparseHandle_t handle = get_handle();
    if (handle == nullptr) return -1;

    PROF_START();

    // --- Upload A (only first call, or if changed) ---
    bool a_was_uploaded = (d_A_row_ptr != nullptr && g_A_n == n && g_A_nnz == nnz_a);
    if (ensure_a_uploaded(n, nnz_a, row_ptr_a, col_ind_a, val_a) != 0)
        return -1;
    if (!a_was_uploaded) {
        PROF_END("upload A (first/changed)");
    }

    // --- Ensure C row_ptr exists (do_spgemm needs it) ---
    if (ensure_c_row_ptr(n) != 0) return -1;

    // --- SpGEMM (algorithm selection + warm-call optimization) ---
    int64_t C_nnz = 0;
    bool success = do_spgemm(handle, n, n, g_desc_A, g_desc_A,
                             &g_spgemm_warm, &C_nnz);

    if (!success) return -1;

#ifdef CUSPARSE_PROFILE
    fprintf(stderr, "[prof r%s]   C_nnz = %lld (full product)\n",
            prof_rank(), (long long)C_nnz);
#endif

    // --- GPU threshold pipeline ---
    return do_threshold_pipeline(n, C_nnz, threshold,
        nnz_c, row_ptr_c_out, col_ind_c_out, val_c_out);
}

// Multi-GPU path: C = A_local * B, then threshold.
// A_local is m × n (row-block of A), B is n × k (column-block of A).
// Both uploaded persistently; SpGEMM + threshold via shared helpers.
int cusparse_gemm_ab_threshold_(
    const int m,            // rows of A_local
    const int n,            // shared dim (cols of A = rows of B)
    const int k,            // cols of B = cols of C
    const int nnz_a,        // nnz of A_local
    const int *row_ptr_a,   // 0-based (NTPoly outer_index)
    const int *col_ind_a,   // 1-based (NTPoly inner_index)
    const double *val_a,
    const int nnz_b,        // nnz of B
    const int *row_ptr_b,   // 0-based
    const int *col_ind_b,   // 1-based
    const double *val_b,
    const double threshold,
    int *nnz_c,
    int **row_ptr_c_out, int **col_ind_c_out, double **val_c_out)
{
    cusparseHandle_t handle = get_handle();
    if (handle == nullptr) return -1;

    PROF_START();

    // --- Upload A_local (first call, or if changed) ---
    bool aloc_was_uploaded = (d_Aloc_row_ptr != nullptr &&
                              g_Aloc_m == m && g_Aloc_nnz == nnz_a);
    if (ensure_aloc_uploaded(m, n, nnz_a, row_ptr_a, col_ind_a, val_a) != 0)
        return -1;
    if (!aloc_was_uploaded)
        PROF_END("upload A_local (first/changed)");

    // --- Upload B (first call, or if changed) ---
    bool b_was_uploaded = (d_B_row_ptr != nullptr &&
                           g_B_n == n && g_B_k == k && g_B_nnz == nnz_b);
    if (ensure_b_uploaded(n, k, nnz_b, row_ptr_b, col_ind_b, val_b) != 0)
        return -1;
    if (!b_was_uploaded)
        PROF_END("upload B (first/changed)");

    // --- Ensure C row_ptr exists (do_spgemm needs it) ---
    if (ensure_c_row_ptr(m) != 0) return -1;

    // --- SpGEMM: C = A_local * B (m × k) ---
    int64_t C_nnz = 0;
    bool success = do_spgemm(handle, m, k, g_desc_Aloc, g_desc_B,
                             &g_spgemm_ab_warm, &C_nnz);

    if (!success) { *nnz_c = 0; return -1; }

#ifdef CUSPARSE_PROFILE
    fprintf(stderr, "[prof r%s]   C_nnz = %lld (full product, A_local*B)\n",
            prof_rank(), (long long)C_nnz);
#endif

    // --- GPU threshold pipeline ---
    return do_threshold_pipeline(m, C_nnz, threshold,
        nnz_c, row_ptr_c_out, col_ind_c_out, val_c_out);
}

// 2.5D k-split path: partial product P = A_local * B (m x k), NO threshold.
// The full (unthresholded) partial is downloaded to persistent host buffers
// so the Fortran side can MPI-exchange it with the other slices; the merge
// + threshold happens afterwards in cusparse_ksplit_merge_.
int cusparse_gemm_ab_partial_(
    const int m,            // rows of A_local
    const int n,            // shared dim (cols of A = rows of B)
    const int k,            // cols of B = cols of C
    const int nnz_a,
    const int *row_ptr_a,   // 0-based (NTPoly outer_index)
    const int *col_ind_a,   // 1-based (NTPoly inner_index)
    const double *val_a,
    const int nnz_b,
    const int *row_ptr_b,   // 0-based
    const int *col_ind_b,   // 1-based
    const double *val_b,
    int *nnz_p,
    int **row_ptr_p_out, int **col_ind_p_out, double **val_p_out)
{
    *nnz_p = 0;
    *row_ptr_p_out = nullptr;
    *col_ind_p_out = nullptr;
    *val_p_out = nullptr;
    if (m <= 0 || n <= 0 || k <= 0 || nnz_a < 0 || nnz_b < 0) return -1;
    if (g_solver_accumulating && (m != g_solver_rows || k != g_solver_cols)) return -1;
    if (nnz_a == 0 || nnz_b == 0) return 0;
    cusparseHandle_t handle = get_handle();
    if (handle == nullptr) { *nnz_p = 0; return -1; }

    PROF_START();

    if (ensure_aloc_uploaded(m, n, nnz_a, row_ptr_a, col_ind_a, val_a) != 0) {
        *nnz_p = 0; return -1;
    }
    if (ensure_b_uploaded(n, k, nnz_b, row_ptr_b, col_ind_b, val_b) != 0) {
        *nnz_p = 0; return -1;
    }
    if (ensure_c_row_ptr(m) != 0) { *nnz_p = 0; return -1; }

    int64_t C_nnz = 0;
    // Changing solver matrices require fresh workspace size queries.
    g_spgemm_ab_warm = false;
    bool success = run_spgemm_general(handle, m, k, g_desc_Aloc, g_desc_B,
        CUSPARSE_SPGEMM_ALG1, &C_nnz, &g_spgemm_ab_warm, true);
    if (!success) return -1;
    if (C_nnz == 0) return 0;

    if (C_nnz > INT_MAX) return -1;
    if (g_solver_accumulating) {
        // Same default stream as SpGEMM: consume C before the next multiply
        // overwrites it. Thresholding occurs only once, after all chunks.
        solver_accumulate_kernel<<<m, 128>>>(m, k, d_C_row_ptr, d_C_col_ind,
                                         d_C_val, d_ks_dense);
        CU_CHECK(cudaGetLastError());
        *nnz_p = (int)C_nnz; // Count only; output pointers remain null.
        return 0;
    }

    // Download the full partial product to persistent host buffers.
    if (h_p_row_n != m) {
        if (h_p_row_ptr) free(h_p_row_ptr);
        h_p_row_ptr = (int *)malloc(sizeof(int) * (size_t)(m + 1));
        h_p_row_n = m;
    }
    if (h_p_cap < (int)C_nnz) {
        int newcap = (int)C_nnz + (int)(C_nnz / 8) + 1024;
        if (h_p_col_ind) free(h_p_col_ind);
        if (h_p_val) free(h_p_val);
        h_p_col_ind = (int *)malloc(sizeof(int) * (size_t)newcap);
        h_p_val = (double *)malloc(sizeof(double) * (size_t)newcap);
        h_p_cap = newcap;
    }
    CU_CHECK(cudaMemcpy(h_p_row_ptr, d_C_row_ptr,
                        sizeof(int) * (size_t)(m + 1), cudaMemcpyDeviceToHost));
    CU_CHECK(cudaMemcpy(h_p_col_ind, d_C_col_ind,
                        sizeof(int) * (size_t)C_nnz, cudaMemcpyDeviceToHost));
    CU_CHECK(cudaMemcpy(h_p_val, d_C_val,
                        sizeof(double) * (size_t)C_nnz, cudaMemcpyDeviceToHost));

    PROF_END("ksplit partial D2H");

    *nnz_p = (int)C_nnz;
    *row_ptr_p_out = h_p_row_ptr;
    *col_ind_p_out = h_p_col_ind;
    *val_p_out = h_p_val;
    return 0;
}

// 2.5D k-split path: merge ps partial products and threshold.
// Inputs: ps CSR partials concatenated on the host.  row_ptr_all holds
// ps*(m+1) entries (slice s at offset s*(m+1), already offset-fixed so
// col_ind_all/val_all are indexed globally).  Output: thresholded merged
// block (m x kk CSR, col_ind 1-based) in the same persistent host buffers
// the threshold pipeline uses.
int cusparse_ksplit_merge_(
    const int ps, const int m, const int kk, const int total_nnz,
    const int *row_ptr_all, const int *col_ind_all, const double *val_all,
    const double threshold, int *nnz_c,
    int **row_ptr_c_out, int **col_ind_c_out, double **val_c_out)
{
    *nnz_c = 0;
    *row_ptr_c_out = nullptr; *col_ind_c_out = nullptr; *val_c_out = nullptr;
    if (ps < 1 || m < 1 || kk < 1 || total_nnz < 0 || threshold < 0) return -1;
    cusparseHandle_t handle = get_handle();
    if (handle == nullptr) return -1;

    PROF_START();

    long long dense_elems = (long long)m * (long long)kk;
    if (dense_elems <= 0 || dense_elems > (long long)12e9) return -1; // ~96 GB guard

    if (g_solver_accumulating) {
        if (m != g_solver_rows || kk != g_solver_cols || ps < 1 || total_nnz != 0) return -1;
        g_solver_accumulating = false;
    } else {
    g_solver_bitmap_valid=false;

    // Upload the concatenated partials.
    int rp_elems = ps * (m + 1);
    if (g_ks_rp_cap < rp_elems) {
        if (d_ks_rp) { cudaFree(d_ks_rp); d_ks_rp = nullptr; }
        CU_CHECK(cudaMalloc(&d_ks_rp, sizeof(int) * (size_t)rp_elems));
        g_ks_rp_cap = rp_elems;
    }
    CU_CHECK(cudaMemcpy(d_ks_rp, row_ptr_all, sizeof(int) * (size_t)rp_elems,
                        cudaMemcpyHostToDevice));
    if (g_ks_cap < total_nnz) {
        int64_t newcap = (int64_t)total_nnz + total_nnz / 8 + 1024;
        if (d_ks_ci) { cudaFree(d_ks_ci); d_ks_ci = nullptr; }
        if (d_ks_v) { cudaFree(d_ks_v); d_ks_v = nullptr; }
        CU_CHECK(cudaMalloc(&d_ks_ci, sizeof(int) * (size_t)newcap));
        CU_CHECK(cudaMalloc(&d_ks_v, sizeof(double) * (size_t)newcap));
        g_ks_cap = newcap;
    }
    if (total_nnz > 0) {
        CU_CHECK(cudaMemcpy(d_ks_ci, col_ind_all, sizeof(int) * (size_t)total_nnz,
                            cudaMemcpyHostToDevice));
        CU_CHECK(cudaMemcpy(d_ks_v, val_all, sizeof(double) * (size_t)total_nnz,
                            cudaMemcpyHostToDevice));
    }
    PROF_END("ksplit merge H2D upload");

    if (g_ks_dense_n != dense_elems) {
        if (d_ks_dense) { cudaFree(d_ks_dense); d_ks_dense = nullptr; }
        CU_CHECK(cudaMalloc(&d_ks_dense, sizeof(double) * (size_t)dense_elems));
        g_ks_dense_n = dense_elems;
    }
    CU_CHECK(cudaMemset(d_ks_dense, 0, sizeof(double) * (size_t)dense_elems));

    // Scatter-add all partials into the dense accumulator.
    dim3 grid(m, ps);
    ksplit_scatter_kernel<<<grid, 128>>>(m, kk, d_ks_rp, d_ks_ci, d_ks_v,
                                         d_ks_dense);
    CU_CHECK(cudaGetLastError());

    } // Legacy host partial upload/scatter; streaming solver is already summed.

    // Per-row counts of entries above threshold.
    if (ensure_row_count(m) != 0) return -1;
    if (ensure_out_row_ptr(m) != 0) return -1;
    ksplit_count_kernel<<<m, 256>>>(d_ks_dense, m, kk, threshold, d_row_count);
    CU_CHECK(cudaGetLastError());
    CU_CHECK(cudaDeviceSynchronize());

    PROF_END("ksplit merge scatter+count");

    // CPU prefix scan over row counts (small: m ints).
    std::vector<int> h_row_count(m);
    CU_CHECK(cudaMemcpy(h_row_count.data(), d_row_count,
                        sizeof(int) * (size_t)m, cudaMemcpyDeviceToHost));
    if (h_out_n < m + 1) {
        if (h_out_row_ptr) free(h_out_row_ptr);
        h_out_row_ptr = (int *)malloc(sizeof(int) * (size_t)(m + 1));
        h_out_n = m + 1;
    }
    h_out_row_ptr[0] = 0;
    int total_keep = 0;
    for (int r = 0; r < m; r++) {
        if (h_row_count[r] > INT_MAX-total_keep) return -1;
        total_keep += h_row_count[r];
        h_out_row_ptr[r + 1] = total_keep;
    }
    *nnz_c = total_keep;

    CU_CHECK(cudaMemcpy(d_out_row_ptr, h_out_row_ptr,
                        sizeof(int) * (size_t)(m + 1), cudaMemcpyHostToDevice));
    if (ensure_out_data(total_keep) != 0) return -1;

    // Compact the dense buffer back to CSR.
    ksplit_compact_kernel<<<m, 32>>>(d_ks_dense, m, kk, threshold,
                                    d_out_row_ptr, d_out_col_ind, d_out_val);
    CU_CHECK(cudaGetLastError());
    CU_CHECK(cudaDeviceSynchronize());

    PROF_END("ksplit merge compact");

    // Download the thresholded merged block (persistent host buffers).
    if (h_out_capacity_nnz < total_keep) {
        if (h_out_col_ind) free(h_out_col_ind);
        if (h_out_val) free(h_out_val);
        h_out_col_ind = (int *)malloc(sizeof(int) * (size_t)total_keep);
        h_out_val = (double *)malloc(sizeof(double) * (size_t)total_keep);
        h_out_capacity_nnz = total_keep;
    }
    CU_CHECK(cudaMemcpy(h_out_col_ind, d_out_col_ind,
                        sizeof(int) * (size_t)total_keep, cudaMemcpyDeviceToHost));
    CU_CHECK(cudaMemcpy(h_out_val, d_out_val,
                        sizeof(double) * (size_t)total_keep, cudaMemcpyDeviceToHost));
    *row_ptr_c_out = h_out_row_ptr;
    *col_ind_c_out = h_out_col_ind;
    *val_c_out = h_out_val;

    PROF_END("ksplit merge D2H");
    return 0;
}

// Experimental dense-accumulator backend: no symbolic product, chunk packing,
// raw partial arrays, or partial-product downloads. Threshold exactly once.
int cusparse_gemm_ab_dense_(const int m, const int n, const int k,
    const int nnz_a, const int *arp, const int *aci, const double *av,
    const int nnz_b, const int *brp, const int *bci, const double *bv,
    const double threshold, int *nnz_c, int **crp, int **cci, double **cv)
{
    *nnz_c = 0; *crp = nullptr; *cci = nullptr; *cv = nullptr;
    if (m < 1 || n < 1 || k < 1 || nnz_a < 0 || nnz_b < 0 || threshold < 0) return -1;
    if (get_handle() == nullptr) return -1;
    if (cusparse_solver_begin_(m,k) != 0) return -1;
    if (nnz_a && nnz_b) {
        if (ensure_aloc_uploaded(m,n,nnz_a,arp,aci,av) != 0 ||
            ensure_b_uploaded(n,k,nnz_b,brp,bci,bv) != 0) return -1;
        solver_dense_product_kernel<<<m,128>>>(k,d_Aloc_row_ptr,d_Aloc_col_ind,
            d_Aloc_val,d_B_row_ptr,d_B_col_ind,d_B_val,d_ks_dense);
        CU_CHECK(cudaGetLastError());
    }
    const int dummy = 0; const double zero = 0;
    return cusparse_ksplit_merge_(1,m,k,0,&dummy,&dummy,&zero,threshold,
                                  nnz_c,crp,cci,cv);
}

void gpu_ops_alloc_(int*,int*,int*,int**,int**,double**,int*);

// Device panel path: inputs use fully device-resident, zero-based CSR.
// Produce the thresholded result directly in an operation slot. Only the
// resulting nnz count crosses to the host; no CSR payload is downloaded.
int cusparse_gemm_device_dense_(int m,int n,int k,int na,const int *arp,
    const int *aci,const double *av,int nb,const int *brp,const int *bci,
    const double *bv,double threshold,int output_slot)
{
    if(m<1 || n<1 || k<1 || na<0 || nb<0 || threshold<0) return -1;
    if(get_handle()==nullptr) return -1;
    const char *option=std::getenv("YUKAWA_BITMAP_ACCUM");
    bool bitmap=option && std::strcmp(option,"1")==0;
    int words=(k-1)/32+1;
    if(bitmap) {
        if(g_solver_accumulating) return -1;
        const long long elements=(long long)m*k;
        if(elements>(long long)12e9) return -1;
        bool reusable=g_solver_bitmap_valid && g_ks_dense_n==elements &&
                      g_solver_bitmap_rows==m && g_solver_bitmap_columns==k;
        if(g_ks_dense_n!=elements) {
            if(d_ks_dense) { CU_CHECK(cudaFree(d_ks_dense)); d_ks_dense=nullptr; g_ks_dense_n=-1; }
            CU_CHECK(cudaMalloc(&d_ks_dense,sizeof(double)*size_t(elements)));
            g_ks_dense_n=elements;
        }
        size_t bytes=size_t(m)*words*sizeof(unsigned);
        if(g_solver_bitmap_capacity<bytes) {
            if(d_solver_bitmap) { CU_CHECK(cudaFree(d_solver_bitmap)); d_solver_bitmap=nullptr; g_solver_bitmap_capacity=0; }
            CU_CHECK(cudaMalloc(&d_solver_bitmap,bytes)); g_solver_bitmap_capacity=bytes;
            reusable=false;
        }
        if(reusable) {
            solver_bitmap_clear_kernel<<<m,128>>>(k,words,d_ks_dense,d_solver_bitmap);
            CU_CHECK(cudaGetLastError());
        } else {
            CU_CHECK(cudaMemset(d_ks_dense,0,sizeof(double)*size_t(elements)));
            CU_CHECK(cudaMemset(d_solver_bitmap,0,bytes));
        }
        g_solver_bitmap_rows=m; g_solver_bitmap_columns=k; g_solver_bitmap_valid=true;
        g_solver_accumulating=true;
    } else if(cusparse_solver_begin_(m,k)!=0) return -1;
    if(na && nb) {
        if(bitmap) solver_bitmap_product_kernel<<<m,128>>>(k,words,arp,aci,av,brp,bci,bv,d_ks_dense,d_solver_bitmap);
        else solver_dense_product_kernel<<<m,128>>>(k,arp,aci,av,brp,bci,bv,d_ks_dense);
        CU_CHECK(cudaGetLastError());
    }
    if(ensure_row_count(m)!=0 || ensure_out_row_ptr(m)!=0) return -1;
    if(bitmap) solver_bitmap_count_kernel<<<m,256>>>(k,words,d_ks_dense,d_solver_bitmap,threshold,d_row_count);
    else ksplit_count_kernel<<<m,256>>>(d_ks_dense,m,k,threshold,d_row_count);
    CU_CHECK(cudaGetLastError());
    int64_t total=thrust::reduce(thrust::device,d_row_count,d_row_count+m,
                                int64_t(0),thrust::plus<int64_t>());
    if(total<0 || total>INT_MAX) return -1;
    thrust::exclusive_scan(thrust::device,d_row_count,d_row_count+m,d_out_row_ptr);
    int nz=int(total),error=-1;
    CU_CHECK(cudaMemcpy(d_out_row_ptr+m,&nz,sizeof(int),cudaMemcpyHostToDevice));
    int *rp=nullptr,*ci=nullptr; double *values=nullptr;
    gpu_ops_alloc_(&output_slot,&m,&nz,&rp,&ci,&values,&error);
    if(error) return -1;
    CU_CHECK(cudaMemcpyAsync(rp,d_out_row_ptr,sizeof(int)*size_t(m+1),cudaMemcpyDeviceToDevice));
    if(nz) {
        if(bitmap) solver_bitmap_compact_kernel<<<m,32>>>(k,words,d_ks_dense,d_solver_bitmap,threshold,d_out_row_ptr,ci,values);
        else ksplit_compact_kernel<<<m,32>>>(d_ks_dense,m,k,threshold,d_out_row_ptr,ci,values);
        CU_CHECK(cudaGetLastError());
    }
    g_solver_accumulating=false;
    return 0;
}

// Memory admission queries do not create a cuSPARSE handle or allocate CSR.
// Device-only queries let the selector reject duplicate rank/GPU mappings
// before initializing an all-rank NCCL communicator.
int gpu_device_info_(int *device,int *count,int8_t *identity) {
    *device=-1; *count=0; std::memset(identity,0,16);
    if(cudaGetDeviceCount(count)!=cudaSuccess || *count<1) return -1;
    const char *local=std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    if(!local) local=std::getenv("SLURM_LOCALID");
    int rank=local?std::atoi(local):0;
    if(rank<0) return -1;
    *device=rank%*count;
    cudaDeviceProp properties;
    if(cudaGetDeviceProperties(&properties,*device)!=cudaSuccess) return -1;
    std::memcpy(identity,properties.uuid.bytes,16); return 0;
}
int gpu_memory_info_(int64_t *free_bytes,int64_t *total_bytes) {
    *free_bytes=0; *total_bytes=0;
    if(!init_gpu()) return -1;
    size_t available=0,total=0;
    if(cudaMemGetInfo(&available,&total)!=cudaSuccess) return -1;
    if(total>size_t(INT64_MAX)) return -1;
    *free_bytes=int64_t(available); *total_bytes=int64_t(total); return 0;
}

// ---------------------------------------------------------------------------
// P1 (GPU_ROADMAP.md): combine ops on CUDA MANAGED memory.
// A "slot" is one local block in CSR living in managed (unified) memory:
//   row_ptr : m+1 entries, 0-based
//   col_ind : nnz entries, 1-based local (1..kk), SORTED within each row,
//             duplicate-free
//   vals    : nnz entries
// Managed memory means CPU and GPU share the same pointer: Fortran fills
// the arrays directly (no upload step), GPU kernels operate in place, and
// any op can fall back to CPU code touching the same arrays. Sync rules:
// cudaDeviceSynchronize() after GPU kernels before CPU access (done inside
// each op), and the CPU->GPU direction relies on UVM page faulting
// (warmup rounds absorb the first-touch cost; prefetch is not required
// for correctness).
// Sortedness/duplicate-freeness is an invariant: the CSC->CSR transpose in
// the mini-app emits ascending columns per row, the merge of two sorted
// rows is sorted, and scale does not reorder. Increment sums duplicates
// on the fly (C + alpha*B); thresholding is left to the caller (v1: 0).
// ---------------------------------------------------------------------------
// Returns a directly usable address without reading the allocation on CPU.
// Full unified-memory systems also support ordinary Fortran allocations.
static const void *import_device_address(const void *pointer,int pageable,cudaError_t &error) {
    if(!pointer) return nullptr;
    cudaPointerAttributes attributes{};
    cudaError_t status=cudaPointerGetAttributes(&attributes,pointer);
    if(status==cudaErrorInvalidValue) {
        (void)cudaGetLastError();
        return pageable?pointer:nullptr;
    }
    if(status!=cudaSuccess) { error=status; return nullptr; }
    if(attributes.type==cudaMemoryTypeManaged) return pointer;
    if(attributes.type==cudaMemoryTypeDevice) {
        int current=-1;
        status=cudaGetDevice(&current);
        if(status!=cudaSuccess) { error=status; return nullptr; }
        if(attributes.device==current) return pointer;
        return nullptr;
    }
    if(attributes.type==cudaMemoryTypeHost) return attributes.devicePointer;
    return pageable?pointer:nullptr;
}

__global__ void import_one_based_offsets(int count,const int *source,int *output) {
    int j=blockIdx.x*blockDim.x+threadIdx.x;
    if(j<count) output[j]=source[j]+1;
}
struct ImportScratch {
    int *offsets=nullptr;
    void *workspace=nullptr;
    ~ImportScratch() { if(workspace) cudaFree(workspace); if(offsets) cudaFree(offsets); }
};

// Synchronous, read-only borrowing of CSC inputs. Output owns its CSR storage.
// 1 declines to the existing host conversion; -1 is an execution error.
int gpu_ops_import_csc_(int *slot,int *rows,int *columns,int *nnz,
    const int *column_ptr,const int *row_indices,const double *values) {
    if(*slot<1 || *rows<1 || *columns<1 || *nnz<0) return -1;
    int device,pageable=0;
    CU_CHECK(cudaGetDevice(&device));
    CU_CHECK(cudaDeviceGetAttribute(&pageable,cudaDevAttrPageableMemoryAccess,device));
    cudaError_t pointer_error=cudaSuccess;
    auto cp=static_cast<const int*>(import_device_address(column_ptr,pageable,pointer_error));
    auto ri=static_cast<const int*>(import_device_address(row_indices,pageable,pointer_error));
    auto input=static_cast<const double*>(import_device_address(values,pageable,pointer_error));
    if(pointer_error!=cudaSuccess) {
        std::fprintf(stderr,"CSC pointer query failed: %s\n",cudaGetErrorString(pointer_error));
        return -1;
    }
    if(!cp || (*nnz && (!ri || !input))) return 1;
    CU_CHECK(cudaDeviceSynchronize());
    int *rp,*ci,error; double *output;
    gpu_ops_alloc_(slot,rows,nnz,&rp,&ci,&output,&error);
    if(error) return -1;
    if(*nnz==0) return 0;
    cusparseHandle_t handle=get_handle();
    if(!handle) return -1;
    ImportScratch scratch;
    size_t offsets_bytes=(size_t(*columns)+1)*sizeof(int);
    CU_CHECK(cudaMalloc(&scratch.offsets,offsets_bytes));
    size_t workspace_bytes=0,free_bytes,total_bytes;
    cusparseStatus_t status=cusparseCsr2cscEx2_bufferSize(handle,*columns,*rows,*nnz,input,
        scratch.offsets,ri,output,rp,ci,CUDA_R_64F,CUSPARSE_ACTION_NUMERIC,CUSPARSE_INDEX_BASE_ONE,
        CUSPARSE_CSR2CSC_ALG1,&workspace_bytes);
    if(status!=CUSPARSE_STATUS_SUCCESS) return -1;
    CU_CHECK(cudaMemGetInfo(&free_bytes,&total_bytes));
    if(workspace_bytes>free_bytes-free_bytes/5) return 1;
    if(workspace_bytes) CU_CHECK(cudaMalloc(&scratch.workspace,workspace_bytes));
    import_one_based_offsets<<<(*columns+256)/256,256>>>(*columns+1,cp,scratch.offsets);
    CU_CHECK(cudaGetLastError());
    status=cusparseCsr2cscEx2(handle,*columns,*rows,*nnz,input,scratch.offsets,ri,output,rp,ci,
        CUDA_R_64F,CUSPARSE_ACTION_NUMERIC,CUSPARSE_INDEX_BASE_ONE,CUSPARSE_CSR2CSC_ALG1,scratch.workspace);
    if(status!=CUSPARSE_STATUS_SUCCESS) return -1;
    subtract_one_kernel<<<(*rows+256)/256,256>>>(rp,*rows+1);
    CU_CHECK(cudaGetLastError());
    CU_CHECK(cudaDeviceSynchronize());
    return 0;
}

#define OPS_SLOTS 4
typedef struct {
    int *rp;
    int *ci;
    double *val;
    int m;
    int nnz;
    int64_t cap_rp;
    int64_t cap_nnz;
} ops_slot_t;
static std::deque<ops_slot_t> g_ops(OPS_SLOTS);
static std::vector<int> g_free_ops;
static int g_owned_live=0;
static int g_ops_dev = -1;  // bound by gpu_ops_init_ (per local rank)

// managed scratch for increment
static int *g_inc_rp = nullptr;      // (m+1)
static int64_t g_inc_rp_cap = 0;
static int *g_inc_ci = nullptr;      // total
static double *g_inc_val = nullptr;  // total
static int64_t g_inc_cap = 0;
static int *g_inc_counts_tmp = nullptr;  // (m) managed, per-increment call

__global__ void ops_scale_kernel(double *val, int nnz, double alpha) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nnz) val[idx] *= alpha;
}

// Count or write a thresholded two-pointer merge of sorted CSR rows.
__global__ void ops_inc_kernel(
    int m, const int *rp_c, const int *ci_c, const double *val_c,
    const int *rp_b, const int *ci_b, const double *val_b,
    double alpha, double threshold, int *counts,
    const int *rp_out, int *ci_out, double *val_out)
{
    int r = blockIdx.x;
    if (r >= m) return;
    int i = rp_c[r], ie = rp_c[r+1], j = rp_b[r], je = rp_b[r+1];
    int out = rp_out ? rp_out[r] : 0;
    while (i < ie || j < je) {
        int col;
        double v;
        if (j == je || (i < ie && ci_c[i] < ci_b[j])) {
            col = ci_c[i]; v = val_c[i++];
        } else if (i == ie || ci_b[j] < ci_c[i]) {
            col = ci_b[j]; v = alpha*val_b[j++];
        } else {
            col = ci_c[i]; v = val_c[i++] + alpha*val_b[j++];
        }
        if (fabs(v) > threshold) {
            if (rp_out) { ci_out[out] = col; val_out[out] = v; }
            ++out;
        }
    }
    if (!rp_out) counts[r] = out;
}

static void ops_check(cudaError_t status) {
    if (status != cudaSuccess) {
        fprintf(stderr, "GPU operation failed: %s\n", cudaGetErrorString(status));
        std::abort();
    }
}

// Gershgorin column contributions, matching NTPoly's column-wise CPU bounds.
// Only O(local columns) host results leave the GPU for the MPI reductions.
__global__ void ops_bounds_kernel(int rows,const int *rp,const int *ci,
    const double *values,int diagonal,double *lower,double *upper) {
    int row=blockIdx.x;
    if(row>=rows) return;
    for(int j=rp[row]+threadIdx.x;j<rp[row+1];j+=blockDim.x) {
        int column=ci[j]-1;
        double value=values[j];
        atomicAdd(lower+column,column==row+diagonal?value:-fabs(value));
        atomicAdd(upper+column,column==row+diagonal?value:fabs(value));
    }
}

int gpu_ops_bounds_(int *slot,int *columns,int *diagonal,double *lower,double *upper) {
    int index=*slot-1;
    if(index<0 || index>=static_cast<int>(g_ops.size()) || *columns<1) return -1;
    auto &matrix=g_ops[index];
    double *device=nullptr;
    size_t bytes=size_t(*columns)*sizeof(double);
    ops_check(cudaMalloc(&device,2*bytes));
    ops_check(cudaMemset(device,0,2*bytes));
    if(matrix.m>0) {
        ops_bounds_kernel<<<matrix.m,128>>>(matrix.m,matrix.rp,matrix.ci,matrix.val,*diagonal,device,device+*columns);
        ops_check(cudaGetLastError());
    }
    ops_check(cudaDeviceSynchronize());
    ops_check(cudaMemcpy(lower,device,bytes,cudaMemcpyDeviceToHost));
    ops_check(cudaMemcpy(upper,device+*columns,bytes,cudaMemcpyDeviceToHost));
    ops_check(cudaFree(device));
    return 0;
}

// Bind this process to its GPU (per local rank). Must be called before
// any managed allocation. Returns the bound device index so the Fortran
// side can point the OpenACC runtime (NTPoly kernels) at the same device.
void gpu_ops_init_(int *err, int *dev_out)
{
    *err = -1;
    *dev_out = -1;
    cusparseHandle_t handle = get_handle();  // runs init_gpu(): cudaSetDevice
    if (handle == nullptr) return;
    if (cudaGetDevice(&g_ops_dev) != cudaSuccess) return;
    *dev_out = g_ops_dev;
    *err = 0;
}

// Library-owned matrices acquire independent slots beyond the legacy four.
// deque keeps existing references stable when another matrix is created.
int gpu_ops_live_count_() { return g_owned_live; }

void gpu_ops_create_(int *slot)
{
    ++g_owned_live;
    if (g_free_ops.empty()) {
        g_ops.push_back(ops_slot_t{});
        *slot = static_cast<int>(g_ops.size());
    } else {
        *slot = g_free_ops.back(); g_free_ops.pop_back();
    }
}
void gpu_ops_destroy_(int *slot)
{
    int s = *slot - 1;
    if (s < OPS_SLOTS || s >= static_cast<int>(g_ops.size())) return;
    ops_check(cudaDeviceSynchronize());
    auto &matrix = g_ops[s];
    if (matrix.rp) ops_check(cudaFree(matrix.rp));
    if (matrix.ci) ops_check(cudaFree(matrix.ci));
    if (matrix.val) ops_check(cudaFree(matrix.val));
    matrix = ops_slot_t{};
    g_free_ops.push_back(*slot); *slot = 0; --g_owned_live;
}

// Allocate (or grow) slot s and hand the managed pointers back to Fortran,
// which fills them directly (no upload step).
void gpu_ops_alloc_(int *slot, int *m, int *nnz,
                    int **rp_out, int **ci_out, double **val_out, int *err)
{
    int s = *slot - 1;
    *err = -1;
    *rp_out = nullptr; *ci_out = nullptr; *val_out = nullptr;
    if (s < 0 || s >= static_cast<int>(g_ops.size()) || *m < 0 || *nnz < 0) return;
    int64_t nnz64 = (*nnz > 0) ? *nnz : 1;
    if (g_ops[s].cap_rp < *m + 1) {
        if (g_ops[s].rp) cudaFree(g_ops[s].rp);
        if (cudaMallocManaged(&g_ops[s].rp, sizeof(int) * (size_t)(*m + 1))
            != cudaSuccess) return;
        g_ops[s].cap_rp = *m + 1;
    }
    if (g_ops[s].cap_nnz < nnz64) {
        if (g_ops[s].ci) cudaFree(g_ops[s].ci);
        if (g_ops[s].val) cudaFree(g_ops[s].val);
        if (cudaMallocManaged(&g_ops[s].ci, sizeof(int) * (size_t)nnz64)
            != cudaSuccess) return;
        if (cudaMallocManaged(&g_ops[s].val, sizeof(double) * (size_t)nnz64)
            != cudaSuccess) return;
        g_ops[s].cap_nnz = nnz64;
    }
    if (*nnz == 0) {
        if (cudaMemset(g_ops[s].rp, 0, sizeof(int)*(size_t)(*m+1)) != cudaSuccess) return;
        if (cudaDeviceSynchronize() != cudaSuccess) return;
    }
    g_ops[s].m = *m;
    g_ops[s].nnz = *nnz;
    *rp_out = g_ops[s].rp;
    *ci_out = g_ops[s].ci;
    *val_out = g_ops[s].val;
    *err = 0;
}

// Scale: vals *= alpha (in place, device kernel).
void gpu_ops_scale_(int *slot, double *alpha)
{
    int s = *slot - 1;
    if (s < 0 || s >= static_cast<int>(g_ops.size()) || g_ops[s].nnz == 0) return;
    int t = 256, b = (g_ops[s].nnz + t - 1) / t;
    ops_scale_kernel<<<b, t>>>(g_ops[s].val, g_ops[s].nnz, *alpha);
    ops_check(cudaGetLastError());
    ops_check(cudaDeviceSynchronize());
}

// Grand sum of the slot's values (scalar only crosses to the host).
void gpu_ops_grandsum_(int *slot, double *sum_out)
{
    int s = *slot - 1;
    *sum_out = 0.0;
    if (s < 0 || s >= static_cast<int>(g_ops.size()) || g_ops[s].nnz == 0) return;
    double *end = g_ops[s].val + (int64_t)g_ops[s].nnz;
    *sum_out = thrust::reduce(thrust::device, g_ops[s].val, end, 0.0,
                              thrust::plus<double>());
}

// Device-to-device copy of one slot into another (managed memcpy).
void gpu_ops_copy_(int *dst, int *src)
{
    int d = *dst - 1, s = *src - 1;
    if (d < 0 || d >= static_cast<int>(g_ops.size()) || s < 0 || s >= static_cast<int>(g_ops.size())) return;
    int m = g_ops[s].m;
    int nnz = g_ops[s].nnz;
    if (g_ops[d].cap_rp < m + 1) {
        if (g_ops[d].rp) cudaFree(g_ops[d].rp);
        if (cudaMallocManaged(&g_ops[d].rp, sizeof(int) * (size_t)(m + 1))
            != cudaSuccess) { fprintf(stderr, "GPU copy allocation failed\n"); std::abort(); }
        g_ops[d].cap_rp = m + 1;
    }
    if (g_ops[d].cap_nnz < (nnz > 0 ? nnz : 1)) {
        if (g_ops[d].ci) cudaFree(g_ops[d].ci);
        if (g_ops[d].val) cudaFree(g_ops[d].val);
        if (cudaMallocManaged(&g_ops[d].ci, sizeof(int) * (size_t)(nnz > 0 ? nnz : 1))
            != cudaSuccess) { fprintf(stderr, "GPU copy allocation failed\n"); std::abort(); }
        if (cudaMallocManaged(&g_ops[d].val, sizeof(double) * (size_t)(nnz > 0 ? nnz : 1))
            != cudaSuccess) { fprintf(stderr, "GPU copy allocation failed\n"); std::abort(); }
        g_ops[d].cap_nnz = (nnz > 0 ? nnz : 1);
    }
    ops_check(cudaMemcpy(g_ops[d].rp, g_ops[s].rp, sizeof(int) * (size_t)(m + 1),
               cudaMemcpyDeviceToDevice));
    if (nnz > 0) {
        ops_check(cudaMemcpy(g_ops[d].ci, g_ops[s].ci, sizeof(int) * (size_t)nnz,
                   cudaMemcpyDeviceToDevice));
        ops_check(cudaMemcpy(g_ops[d].val, g_ops[s].val, sizeof(double) * (size_t)nnz,
                   cudaMemcpyDeviceToDevice));
    }
    ops_check(cudaDeviceSynchronize());
    g_ops[d].m = m;
    g_ops[d].nnz = nnz;
}

// Increment: C += alpha*B (both slots managed, rows sorted, duplicates
// summed on the fly). Merges into managed scratch, then ping-pongs the
// scratch into slot C. The CPU does the row_ptr prefix scan directly on
// the managed counts (no explicit D2H copy). Threshold after summation.
void gpu_ops_increment_threshold_(int *slot_c, int *slot_b, double *alpha,
                        double *threshold, int *nnz_out)
{
    *nnz_out = -1;
    if (*threshold < 0) return;
    int c = *slot_c - 1, b = *slot_b - 1;
    if (c < 0 || c >= static_cast<int>(g_ops.size()) || b < 0 || b >= static_cast<int>(g_ops.size())) { *nnz_out = -1; return; }
    int m = g_ops[c].m;
    if (m != g_ops[b].m) { *nnz_out = -1; return; }
    if (m == 0) { *nnz_out = 0; return; }

    if (g_inc_rp_cap < m + 1) {
        if (g_inc_rp) cudaFree(g_inc_rp);
        if (cudaMallocManaged(&g_inc_rp, sizeof(int) * (size_t)(m + 1))
            != cudaSuccess) { *nnz_out = -1; return; }
        g_inc_rp_cap = m + 1;
    }
    if (cudaMallocManaged(&g_inc_counts_tmp, sizeof(int) * (size_t)m)
        != cudaSuccess) { *nnz_out = -1; return; }

    // Pass 1: per-row merged counts.
    ops_inc_kernel<<<m, 1>>>(m, g_ops[c].rp, g_ops[c].ci, g_ops[c].val,
        g_ops[b].rp, g_ops[b].ci, g_ops[b].val, *alpha, *threshold,
        g_inc_counts_tmp, nullptr, nullptr, nullptr);
    ops_check(cudaGetLastError());
    ops_check(cudaDeviceSynchronize());

    // CPU prefix scan directly on managed memory.
    g_inc_rp[0] = 0;
    int total = 0;
    for (int r = 0; r < m; r++) {
        if (g_inc_counts_tmp[r] > INT_MAX-total) { *nnz_out = -1; return; }
        total += g_inc_counts_tmp[r];
        g_inc_rp[r + 1] = total;
    }
    cudaFree(g_inc_counts_tmp);

    // Pass 2: merge-write into managed scratch.
    if (g_inc_cap < (total > 0 ? total : 1)) {
        if (g_inc_ci) cudaFree(g_inc_ci);
        if (g_inc_val) cudaFree(g_inc_val);
        if (cudaMallocManaged(&g_inc_ci, sizeof(int) * (size_t)(total > 0 ? total : 1))
            != cudaSuccess) { *nnz_out = -1; return; }
        if (cudaMallocManaged(&g_inc_val, sizeof(double) * (size_t)(total > 0 ? total : 1))
            != cudaSuccess) { *nnz_out = -1; return; }
        g_inc_cap = (total > 0 ? total : 1);
    }
    ops_inc_kernel<<<m, 1>>>(m, g_ops[c].rp, g_ops[c].ci, g_ops[c].val,
        g_ops[b].rp, g_ops[b].ci, g_ops[b].val, *alpha, *threshold,
        nullptr, g_inc_rp, g_inc_ci, g_inc_val);
    ops_check(cudaGetLastError());
    ops_check(cudaDeviceSynchronize());

    // Ping-pong: slot C takes the scratch, scratch keeps the old buffers.
    int *t_rp = g_ops[c].rp; g_ops[c].rp = g_inc_rp; g_inc_rp = t_rp;
    int *t_ci = g_ops[c].ci; g_ops[c].ci = g_inc_ci; g_inc_ci = t_ci;
    double *t_val = g_ops[c].val; g_ops[c].val = g_inc_val; g_inc_val = t_val;
    std::swap(g_ops[c].cap_rp, g_inc_rp_cap);
    std::swap(g_ops[c].cap_nnz, g_inc_cap);
    g_ops[c].m = m;
    g_ops[c].nnz = total;
    *nnz_out = total;
}

// Preserve the P1 interface; zero threshold also removes exact cancellations.
void gpu_ops_increment_(int *slot_c, int *slot_b, double *alpha, int *nnz_out)
{
    double threshold = 0.0;
    gpu_ops_increment_threshold_(slot_c, slot_b, alpha, &threshold, nnz_out);
}

// --- P2 reductions (GPU_ROADMAP.md): trace / dot / column abs-sums ---
// All compute a per-block scalar (or per-column vector) on device; the
// MPI collectives that combine ranks stay on the Fortran side, which
// owns the grid communicators. NTPoly semantics (verified against
// Source/Fortran/distributed_algebra_includes):
//   trace    : sum of entries with global row == global col
//   dot      : sum over common (r,c) of A*B (Hadamard then grand sum)
//   colsum   : per-column sum of |v| (feeds MatrixNorm's max column
//              1-norm and MatrixSigma's 1/max^2)

// One thread per row; rows are sorted so the diagonal entry (if any) is
// found by a short scan. diag_offset = (my_row - my_col) * m_loc: the
// block-local column index that corresponds to the global diagonal
// (off-diagonal blocks get an offset outside [1, kk] and contribute 0).
__global__ void ops_trace_kernel(int m, const int *rp, const int *ci,
                                 const double *val, int diag_offset,
                                 double *row_sums)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= m) return;
    double local = 0.0;
    int target = r + 1 + diag_offset;
    if (target >= 1) {
        int lo = rp[r], hi = rp[r + 1];
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            if (ci[mid] < target) lo = mid + 1;
            else hi = mid;
        }
        if (lo < rp[r + 1] && ci[lo] == target) local = val[lo];
    }
    row_sums[r] = local;
}

// One thread per row; two-pointer walk over the two sorted rows, summing
// products where the column matches.
__global__ void ops_dot_kernel(int m, const int *rp_a, const int *ci_a,
                               const double *va, const int *rp_b,
                               const int *ci_b, const double *vb,
                               double *row_sums)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= m) return;
    int i = rp_a[r], i_end = rp_a[r + 1];
    int j = rp_b[r], j_end = rp_b[r + 1];
    double local = 0.0;
    while (i < i_end && j < j_end) {
        int ca = ci_a[i], cb = ci_b[j];
        if (ca == cb) { local += va[i] * vb[j]; ++i; ++j; }
        else if (ca < cb) ++i;
        else ++j;
    }
    row_sums[r] = local;
}

// Per-column sum of |v| for one slot, written into managed scratch.
// CSR storage means the column accumulation uses atomicAdd; the sums are
// short (sparse columns) so ordering noise is ~1e-16 relative.
__global__ void ops_colsum_kernel(int m, int kk, const int *rp,
                                  const int *ci, const double *val,
                                  double *colsum)
{
    int r = blockIdx.x;
    if (r >= m) return;
    for (int j = rp[r] + threadIdx.x; j < rp[r + 1]; j += blockDim.x)
        atomicAdd(&colsum[ci[j] - 1], fabs(val[j]));
}

// Trace of the slot's block. diag_offset = (my_row - my_col) * m_loc.
void gpu_ops_trace_(int *slot, int *diag_offset, double *trace_out)
{
    int s = *slot - 1;
    *trace_out = 0.0;
    if (s < 0 || s >= static_cast<int>(g_ops.size()) || g_ops[s].nnz == 0) return;
    int m = g_ops[s].m;
    double *d_rows = nullptr;
    if (cudaMallocManaged(&d_rows, sizeof(double) * (size_t)m)
        != cudaSuccess) { fprintf(stderr, "GPU reduction allocation failed\n"); std::abort(); }
    int t = 128, b = (m + t - 1) / t;
    ops_trace_kernel<<<b, t>>>(m, g_ops[s].rp, g_ops[s].ci, g_ops[s].val,
                               *diag_offset, d_rows);
    double *end = d_rows + m;
    *trace_out = thrust::reduce(thrust::device, d_rows, end, 0.0,
                                thrust::plus<double>());
    cudaFree(d_rows);
    ops_check(cudaGetLastError());
    ops_check(cudaDeviceSynchronize());
}

// Dot product of two slots over their common entries.
void gpu_ops_dot_(int *slot_a, int *slot_b, double *dot_out)
{
    int a = *slot_a - 1, b = *slot_b - 1;
    *dot_out = 0.0;
    if (a < 0 || a >= static_cast<int>(g_ops.size()) || b < 0 || b >= static_cast<int>(g_ops.size())) return;
    int m = g_ops[a].m;
    if (m != g_ops[b].m || g_ops[a].nnz == 0 || g_ops[b].nnz == 0) return;
    double *d_rows = nullptr;
    if (cudaMallocManaged(&d_rows, sizeof(double) * (size_t)m)
        != cudaSuccess) { fprintf(stderr, "GPU reduction allocation failed\n"); std::abort(); }
    int t = 128, blocks = (m + t - 1) / t;
    ops_dot_kernel<<<blocks, t>>>(m, g_ops[a].rp, g_ops[a].ci, g_ops[a].val,
                                  g_ops[b].rp, g_ops[b].ci, g_ops[b].val,
                                  d_rows);
    double *end = d_rows + m;
    *dot_out = thrust::reduce(thrust::device, d_rows, end, 0.0,
                              thrust::plus<double>());
    cudaFree(d_rows);
    ops_check(cudaGetLastError());
    ops_check(cudaDeviceSynchronize());
}

// Per-column |v| sums of one slot into a HOST array (kk entries), for the
// Fortran-side column_comm/row_comm MPI collectives (MatrixNorm,
// MatrixSigma).
void gpu_ops_colsum_(int *slot, int *kk, double *colsum_out)
{
    int s = *slot - 1;
    for (int j = 0; j < *kk; ++j) colsum_out[j] = 0.0;
    if (s < 0 || s >= static_cast<int>(g_ops.size()) || g_ops[s].nnz == 0) return;
    int m = g_ops[s].m;
    double *d_cols = nullptr;
    if (cudaMallocManaged(&d_cols, sizeof(double) * (size_t)(*kk))
        != cudaSuccess) { fprintf(stderr, "GPU reduction allocation failed\n"); std::abort(); }
    ops_check(cudaMemset(d_cols, 0, sizeof(double) * (size_t)(*kk)));
    ops_colsum_kernel<<<m, 128>>>(m, *kk, g_ops[s].rp, g_ops[s].ci,
                                  g_ops[s].val, d_cols);
    ops_check(cudaGetLastError());
    ops_check(cudaDeviceSynchronize());
    ops_check(cudaMemcpy(colsum_out, d_cols, sizeof(double) * (size_t)(*kk),
               cudaMemcpyDeviceToHost));
    cudaFree(d_cols);
}

// Hand the slot's managed pointers back to Fortran (the CPU reads them
// directly after a sync — the managed-memory payoff: no download step).
void gpu_ops_get_(int *slot, int **rp, int **ci, double **val,
                  int *m, int *nnz)
{
    int s = *slot - 1;
    if (s < 0 || s >= static_cast<int>(g_ops.size())) { *m = -1; *nnz = -1; return; }
    ops_check(cudaDeviceSynchronize());
    *rp = g_ops[s].rp;
    *ci = g_ops[s].ci;
    *val = g_ops[s].val;
    *m = g_ops[s].m;
    *nnz = g_ops[s].nnz;
}

// Reuse the sorted merge/compaction path with a separate empty source.
// A self increment with alpha=0 would evaluate 0*Inf and change values.
void gpu_ops_prune_(int *slot, double *threshold, int *nnz)
{
    *nnz = -1;
    int s = *slot - 1;
    if (s < 0 || s >= static_cast<int>(g_ops.size()) || !std::isfinite(*threshold) || *threshold < 0) return;
    int rows = g_ops[s].m, empty = 0, zero = 0, status = -1;
    int *rp = nullptr, *ci = nullptr;
    double *values = nullptr, alpha = 0;
    gpu_ops_create_(&empty);
    gpu_ops_alloc_(&empty, &rows, &zero, &rp, &ci, &values, &status);
    if (status == 0) gpu_ops_increment_threshold_(slot, &empty, &alpha, threshold, nnz);
    gpu_ops_destroy_(&empty);
}

void cusparse_gemm_free_(int *row_ptr, int *col_ind, double *val) {
    // No-op: buffers are persistent (freed at program exit).
}

} // extern "C"
