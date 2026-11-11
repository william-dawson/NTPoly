// Minimal descriptor API for exercising the real upload helpers on the host.
constexpr int cudaMemcpyHostToDevice=0;
struct Descriptor { int rows,columns,nnz; };
using cusparseSpMatDescr_t=Descriptor*;
constexpr int CUSPARSE_INDEX_32I=0, CUSPARSE_INDEX_BASE_ZERO=0, CUDA_R_64F=0;
inline int cusparseDestroySpMat(Descriptor *p) { delete p; return 0; }
inline int cusparseCreateCsr(Descriptor **p,int m,int n,int nz,int*,int*,double*,int,int,int,int) {
    *p=new Descriptor{m,n,nz}; return 0;
}
#define CU_CHECK(call) do { if((call)!=0) return -1; } while(0)
#define CS_CHECK(call) CU_CHECK(call)
