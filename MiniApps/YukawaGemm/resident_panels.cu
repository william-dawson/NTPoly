// GPU CSR panel construction for the one-node purification experiment.
// Local slots: zero-based row pointers, one-based sorted column indices.
// Assembled multiply panels: entirely device allocated, zero-based CSR.
#include <cuda_runtime.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdio>
#include <vector>

extern "C" {
void gpu_ops_get_(int*,int**,int**,double**,int*,int*);
int solver_nccl_group_info_(int,int*,int*);
int solver_nccl_gather_i_(const int*,const int*,const int*,int*,const int*,const int*);
int solver_nccl_device_gather_(int,const void*,int,void*,const int*,const int*,int);
int cusparse_gemm_device_dense_(int,int,int,int,const int*,const int*,const double*,
                               int,const int*,const int*,const double*,double,int);
}
#define DCUDA(call) do { auto e=(call); if(e!=cudaSuccess) { \
    std::fprintf(stderr,"Resident panels %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e));return -1;} }while(0)
#define CHECK(call) do { if((call)!=0) return -1; }while(0)
namespace {
struct Buffer {
    void *data=nullptr; size_t capacity=0;
    int reserve(size_t bytes) {
        bytes=std::max(bytes,size_t(8));
        if(bytes<=capacity) return 0;
        if(data) { DCUDA(cudaFree(data)); data=nullptr; capacity=0; }
        DCUDA(cudaMalloc(&data,bytes)); capacity=bytes; return 0;
    }
    template<class T> T *as() { return static_cast<T*>(data); }
    int clear() { if(data) DCUDA(cudaFree(data)); data=nullptr; capacity=0; return 0; }
};
struct CSR { Buffer rp,ci,values; } local,row_received,column_received,a,b;
Buffer row_displacements,column_displacements;
double elapsed(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
__global__ void horizontal_counts(int rows,int peers,const int *peer_rp,int *rp) {
    int row=blockIdx.x*blockDim.x+threadIdx.x;
    if(row==0) rp[0]=0;
    if(row>=rows) return;
    int count=0;
    for(int peer=0;peer<peers;++peer) {
        const int *p=peer_rp+size_t(peer)*(rows+1);
        count+=p[row+1]-p[row];
    }
    rp[row+1]=count;
}
__global__ void horizontal_pack(int rows,int columns,int peers,const int *peer_rp,
    const int *peer_ci,const double *peer_values,const int *displs,
    const int *rp,int *ci,double *values) {
    int row=blockIdx.x,output=rp[row];
    for(int peer=0;peer<peers;++peer) {
        const int *p=peer_rp+size_t(peer)*(rows+1);
        for(int j=p[row]+threadIdx.x;j<p[row+1];j+=blockDim.x) {
            int source=displs[peer]+j,destination=output+j-p[row];
            ci[destination]=peer*columns+peer_ci[source]-1;
            values[destination]=peer_values[source];
        }
        output+=p[row+1]-p[row];
    }
}
__global__ void vertical_rows(int rows,int peers,const int *peer_rp,
                              const int *displs,int total,int *rp) {
    int global=blockIdx.x*blockDim.x+threadIdx.x;
    if(global==rows*peers) rp[global]=total;
    else if(global<rows*peers) {
        int peer=global/rows,row=global%rows;
        rp[global]=displs[peer]+peer_rp[size_t(peer)*(rows+1)+row];
    }
}
__global__ void zero_base_columns(int count,int *ci) {
    int j=blockIdx.x*blockDim.x+threadIdx.x;
    if(j<count) --ci[j];
}
int exchange(int which,int rows,int nnz,CSR &received,Buffer &device_displs,
             std::vector<int> &counts,std::vector<int> &displs,int &total) {
    int peers,rank;
    CHECK(solver_nccl_group_info_(which,&peers,&rank));
    counts.resize(peers); displs.resize(peers);
    std::vector<int> ones(peers,1),offsets(peers);
    for(int p=0;p<peers;++p) offsets[p]=p;
    int one=1;
    CHECK(solver_nccl_gather_i_(&which,&nnz,&one,counts.data(),ones.data(),offsets.data()));
    total=0;
    for(int p=0;p<peers;++p) {
        if(counts[p]<0 || counts[p]>INT_MAX-total) return -1;
        displs[p]=total; total+=counts[p];
    }
    CHECK(received.rp.reserve(size_t(peers)*(rows+1)*sizeof(int)));
    CHECK(received.ci.reserve(size_t(total)*sizeof(int)));
    CHECK(received.values.reserve(size_t(total)*sizeof(double)));
    CHECK(device_displs.reserve(size_t(peers)*sizeof(int)));
    DCUDA(cudaMemcpyAsync(device_displs.data,displs.data(),size_t(peers)*sizeof(int),cudaMemcpyHostToDevice));
    std::vector<int> row_counts(peers,rows+1),row_displs(peers);
    for(int p=0;p<peers;++p) row_displs[p]=p*(rows+1);
    CHECK(solver_nccl_device_gather_(which,local.rp.data,rows+1,received.rp.data,
                                    row_counts.data(),row_displs.data(),0));
    CHECK(solver_nccl_device_gather_(which,local.ci.data,nnz,received.ci.data,
                                    counts.data(),displs.data(),0));
    CHECK(solver_nccl_device_gather_(which,local.values.data,nnz,received.values.data,
                                    counts.data(),displs.data(),1));
    // Keep pageable metadata alive until its small H2D copy has completed.
    DCUDA(cudaDeviceSynchronize());
    return 0;
}
int snapshot(int slot,int rows,int &nnz) {
    int *rp,*ci,slot_rows; double *values;
    gpu_ops_get_(&slot,&rp,&ci,&values,&slot_rows,&nnz);
    if(slot_rows!=rows || nnz<0) return -1;
    CHECK(local.rp.reserve(size_t(rows+1)*sizeof(int)));
    CHECK(local.ci.reserve(size_t(nnz)*sizeof(int)));
    CHECK(local.values.reserve(size_t(nnz)*sizeof(double)));
    DCUDA(cudaMemcpyAsync(local.rp.data,rp,size_t(rows+1)*sizeof(int),cudaMemcpyDeviceToDevice));
    if(nnz) {
        DCUDA(cudaMemcpyAsync(local.ci.data,ci,size_t(nnz)*sizeof(int),cudaMemcpyDeviceToDevice));
        DCUDA(cudaMemcpyAsync(local.values.data,values,size_t(nnz)*sizeof(double),cudaMemcpyDeviceToDevice));
    }
    return 0;
}
}

extern "C" int gpu_resident_multiply_(int left_slot,int right_slot,int output_slot,int rows,int n,
    int columns,int process_rows,int process_columns,double threshold,double *phases) {
    if(rows<1 || columns<1 || n<1 || process_rows<1 || process_columns<1 ||
       left_slot==output_slot || right_slot==output_slot ||
       (long long)rows*process_rows!=n || (long long)columns*process_columns!=n) return -1;
    int peers,rank;
    CHECK(solver_nccl_group_info_(1,&peers,&rank));
    if(peers!=process_columns) return -1;
    CHECK(solver_nccl_group_info_(2,&peers,&rank));
    if(peers!=process_rows) return -1;
    auto start=std::chrono::steady_clock::now();
    int nnz;
    CHECK(snapshot(left_slot,rows,nnz));
    std::vector<int> counts,displs;
    int nnz_a,nnz_b;
    CHECK(exchange(1,rows,nnz,row_received,row_displacements,counts,displs,nnz_a));
    phases[0]=elapsed(start); start=std::chrono::steady_clock::now();
    CHECK(a.rp.reserve(size_t(rows+1)*sizeof(int)));
    CHECK(a.ci.reserve(size_t(nnz_a)*sizeof(int)));
    CHECK(a.values.reserve(size_t(nnz_a)*sizeof(double)));
    horizontal_counts<<<(rows+255)/256,256>>>(rows,process_columns,row_received.rp.as<int>(),a.rp.as<int>());
    DCUDA(cudaGetLastError());
    thrust::inclusive_scan(thrust::device,a.rp.as<int>(),a.rp.as<int>()+rows+1,a.rp.as<int>());
    horizontal_pack<<<rows,128>>>(rows,columns,process_columns,row_received.rp.as<int>(),
        row_received.ci.as<int>(),row_received.values.as<double>(),row_displacements.as<int>(),
        a.rp.as<int>(),a.ci.as<int>(),a.values.as<double>());
    DCUDA(cudaGetLastError()); DCUDA(cudaDeviceSynchronize());
    phases[1]=elapsed(start); start=std::chrono::steady_clock::now();
    if(left_slot!=right_slot) CHECK(snapshot(right_slot,rows,nnz));
    CHECK(exchange(2,rows,nnz,column_received,column_displacements,counts,displs,nnz_b));
    CHECK(b.rp.reserve(size_t(n+1)*sizeof(int)));
    vertical_rows<<<(n+1+255)/256,256>>>(rows,process_rows,column_received.rp.as<int>(),
        column_displacements.as<int>(),nnz_b,b.rp.as<int>());
    if(nnz_b) zero_base_columns<<<(nnz_b+255)/256,256>>>(nnz_b,column_received.ci.as<int>());
    DCUDA(cudaGetLastError()); DCUDA(cudaDeviceSynchronize());
    phases[2]=elapsed(start); start=std::chrono::steady_clock::now();
    CHECK(cusparse_gemm_device_dense_(rows,n,columns,nnz_a,a.rp.as<int>(),a.ci.as<int>(),
        a.values.as<double>(),nnz_b,b.rp.as<int>(),column_received.ci.as<int>(),
        column_received.values.as<double>(),threshold,output_slot));
    DCUDA(cudaDeviceSynchronize()); phases[3]=elapsed(start);
    return 0;
}

extern "C" int gpu_resident_square_(int input_slot,int output_slot,int rows,int n,
    int columns,int process_rows,int process_columns,double threshold,double *phases) {
    return gpu_resident_multiply_(input_slot,input_slot,output_slot,rows,n,columns,
                                  process_rows,process_columns,threshold,phases);
}

extern "C" int gpu_resident_free_() {
    DCUDA(cudaDeviceSynchronize());
    for(CSR *p:{&local,&row_received,&column_received,&a,&b}) {
        CHECK(p->rp.clear()); CHECK(p->ci.clear()); CHECK(p->values.clear());
    }
    CHECK(row_displacements.clear()); CHECK(column_displacements.clear()); return 0;
}
