// Real CUDA test: caller-owned CSC input is built/changed on GPU and borrowed.
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
extern "C" {
void gpu_ops_init_(int*,int*);
void gpu_ops_create_(int*);
void gpu_ops_destroy_(int*);
int gpu_ops_import_csc_(int*,int*,int*,int*,const int*,const int*,const double*);
int gpu_ops_export_csc_(int*,int*,int*,int*,double*);
void gpu_ops_get_(int*,int**,int**,double**,int*,int*);
}
static void check(bool ok) { if(!ok) { std::fprintf(stderr,"CSC import test failed\n"); std::exit(1); } }
static void cuda_check(cudaError_t s) { check(s==cudaSuccess); }
__global__ void change_values(int n,double *values,double factor) {
    int j=blockIdx.x*blockDim.x+threadIdx.x;
    if(j<n) values[j]*=factor;
}
static void *allocate(size_t bytes,int mode) {
    void *pointer=nullptr; bytes=std::max(bytes,size_t(1));
    if(mode==0) cuda_check(cudaMallocManaged(&pointer,bytes));
    if(mode==1) cuda_check(cudaMalloc(&pointer,bytes));
    if(mode==2) { pointer=std::malloc(bytes); check(pointer); }
    if(mode==3) cuda_check(cudaHostAlloc(&pointer,bytes,cudaHostAllocMapped));
    return pointer;
}
static void release(void *pointer,int mode) {
    if(mode==2) std::free(pointer);
    else if(mode==3) cuda_check(cudaFreeHost(pointer));
    else cuda_check(cudaFree(pointer));
}
int main() {
    int error,device,slot;
    gpu_ops_init_(&error,&device); check(error==0);
    gpu_ops_create_(&slot);
    int pageable=0;
    cuda_check(cudaDeviceGetAttribute(&pageable,cudaDevAttrPageableMemoryAccess,device));
    cudaStream_t stream; cuda_check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    int borrowed=0,declined=0;
    for(int mode=0;mode<4;++mode) for(int repeat=0;repeat<6;++repeat) {
        int rows=repeat%2?7:3,columns=repeat%3?5:9;
        std::vector<int> cp(columns+1,0),ri;
        std::vector<double> values,dense(rows*columns,0);
        for(int c=0;c<columns;++c) {
            for(int r=0;r<rows;++r) if(repeat!=2 && (r*7+c*11+repeat)%5==0) {
                double v=(r+1)*.25-(c+1)*.1;
                ri.push_back(r+1); values.push_back(v); dense[r*columns+c]=v;
            }
            cp[c+1]=ri.size();
        }
        int nnz=values.size();
        auto dcp=static_cast<int*>(allocate(cp.size()*sizeof(int),mode));
        auto dri=static_cast<int*>(allocate(ri.size()*sizeof(int),mode));
        auto dv=static_cast<double*>(allocate(values.size()*sizeof(double),mode));
        cuda_check(cudaMemcpy(dcp,cp.data(),cp.size()*sizeof(int),cudaMemcpyHostToDevice));
        if(nnz) {
            cuda_check(cudaMemcpy(dri,ri.data(),ri.size()*sizeof(int),cudaMemcpyHostToDevice));
            cuda_check(cudaMemcpy(dv,values.data(),values.size()*sizeof(double),cudaMemcpyHostToDevice));
        }
        // The producer runs on another stream. Import must observe its work,
        // without a CPU payload read, including for managed inputs.
        bool accessible=mode!=2 || pageable;
        if(accessible && nnz) {
            change_values<<<(nnz+255)/256,256,0,stream>>>(nnz,dv,2.0);
            cuda_check(cudaGetLastError());
        }
        int status=gpu_ops_import_csc_(&slot,&rows,&columns,&nnz,dcp,dri,dv);
        if(!accessible) { check(status==1); ++declined; }
        else {
            check(status==0); ++borrowed;
            int *rp,*ci,m,n; double *output;
            gpu_ops_get_(&slot,&rp,&ci,&output,&m,&n);
            check(m==rows && n==nnz && rp[0]==0 && rp[rows]==nnz);
            for(int r=0;r<rows;++r) {
                int expected=0;
                for(int c=0;c<columns;++c) if(repeat!=2 && (r*7+c*11+repeat)%5==0) {
                    int j=rp[r]+expected++;
                    check(j<rp[r+1] && ci[j]==c+1 && std::abs(output[j]-2*dense[r*columns+c])<1e-14);
                }
                check(expected==rp[r+1]-rp[r]);
            }
            if(nnz) {
                change_values<<<(nnz+255)/256,256,0,stream>>>(nnz,dv,3.0);
                cuda_check(cudaGetLastError());
                check(gpu_ops_import_csc_(&slot,&rows,&columns,&nnz,dcp,dri,dv)==0);
                gpu_ops_get_(&slot,&rp,&ci,&output,&m,&n);
                for(int r=0;r<rows;++r) for(int j=rp[r];j<rp[r+1];++j)
                    check(std::abs(output[j]-6*dense[r*columns+ci[j]-1])<1e-14);
            }
            // Round-trip all destination allocation kinds, including changed
            // source values, empty matrices, stored zeros and rectangular shapes.
            check(gpu_ops_export_csc_(&slot,&columns,nullptr,nullptr,nullptr)==1);
            for(int destination_mode=0;destination_mode<4;++destination_mode) {
                auto ocp=static_cast<int*>(allocate(cp.size()*sizeof(int),destination_mode));
                auto ori=static_cast<int*>(allocate(ri.size()*sizeof(int),destination_mode));
                auto ov=static_cast<double*>(allocate(values.size()*sizeof(double),destination_mode));
                int exported=gpu_ops_export_csc_(&slot,&columns,ocp,ori,ov);
                if(destination_mode==2 && !pageable)check(exported==1);
                else {
                    check(exported==0);
                    std::vector<int> hcp(cp.size()),hri(ri.size());std::vector<double> hv(values.size());
                    cuda_check(cudaMemcpy(hcp.data(),ocp,hcp.size()*sizeof(int),cudaMemcpyDeviceToHost));
                    check(hcp==cp);
                    if(nnz) {
                        cuda_check(cudaMemcpy(hri.data(),ori,hri.size()*sizeof(int),cudaMemcpyDeviceToHost));
                        cuda_check(cudaMemcpy(hv.data(),ov,hv.size()*sizeof(double),cudaMemcpyDeviceToHost));
                        check(hri==ri);
                        for(int j=0;j<nnz;++j)check(std::abs(hv[j]-6*values[j])<1e-14);
                    }
                }
                release(ocp,destination_mode);release(ori,destination_mode);release(ov,destination_mode);
            }
        }
        release(dcp,mode); release(dri,mode); release(dv,mode);
    }
    gpu_ops_destroy_(&slot); check(slot==0);
    cuda_check(cudaStreamDestroy(stream));
    std::printf("Borrowed CSC import/export: managed/device/system/mapped, changed values, "
                "rectangular/empty/reuse/sorted PASS borrowed=%d declined=%d pageable=%d\n",
                borrowed,declined,pageable);
}
