// Test-only serial CUDA API/kernel emulation. This does not validate CUDA.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <climits>
#include <utility>
#include <functional>
#include <cassert>
#include <vector>
#include <deque>
#include <algorithm>
#define __global__
struct Index { int x=0; };
static Index blockIdx, threadIdx, blockDim;
template<class F> void launch(int blocks, int threads, F f) {
    blockDim.x=threads;
    for(blockIdx.x=0; blockIdx.x<blocks; ++blockIdx.x)
        for(threadIdx.x=0; threadIdx.x<threads; ++threadIdx.x) f();
}
using cudaError_t=int;
constexpr int cudaSuccess=0, cudaMemcpyDeviceToDevice=0, cudaMemcpyDeviceToHost=0;
template<class T> int cudaMallocManaged(T **p, size_t n) {
    *p=static_cast<T*>(malloc(n)); return *p ? 0 : 1;
}
inline int cudaFree(void *p) { free(p); return 0; }
inline int cudaMemset(void *p, int v, size_t n) { memset(p,v,n); return 0; }
inline int cudaMemcpy(void *d,const void *s,size_t n,int) { memmove(d,s,n); return 0; }
inline int cudaDeviceSynchronize() { return 0; }
inline int cudaGetLastError() { return 0; }
inline int cudaGetDevice(int *d) { *d=0; return 0; }
inline const char *cudaGetErrorString(int) { return "host test CUDA error"; }
inline void atomicAdd(double *p,double v) { *p+=v; }
using cusparseHandle_t=void*;
inline cusparseHandle_t get_handle() { return reinterpret_cast<void*>(1); }
namespace thrust {
constexpr int device=0;
template<class T> using plus=std::plus<T>;
template<class P,class T,class F> T reduce(int,P p,P end,T sum,F f) {
    for(;p!=end;++p) sum=f(sum,*p);
    return sum;
}
}
