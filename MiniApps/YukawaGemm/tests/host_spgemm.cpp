// Deterministic CPU substitutes for cuSPARSE, solely for driver/MPI tests.
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
static bool accumulating = false;
static std::vector<double> accumulated;
static int accumulated_m = 0, accumulated_k = 0;
extern "C" {
void gpu_profile_begin_() {}
void gpu_profile_end_() {}
int gpu_resident_square_() { std::abort(); }
int gpu_resident_free_() { return 0; }
int cusparse_solver_begin_(int m,int k) {
    if(accumulating || m<1 || k<1) return -1;
    accumulating=true; accumulated_m=m; accumulated_k=k;
    accumulated.assign((size_t)m*k,0); return 0;
}
int cusparse_gemm_ab_partial_(int m,int n,int k,int na,const int *ar,
    const int *ai,const double *av,int nb,const int *br,const int *bi,
    const double *bv,int *nz,int **rp,int **ci,double **val) {
    if(getenv("NTG_TEST_FAIL_SPGEMM")) return -1;
    static std::vector<int> r,c;
    static std::vector<double> v;
    r.assign(m+1,0); c.clear(); v.clear();
    for(int row=0;row<m;++row) {
        std::vector<double> dense(k,0);
        for(int a=ar[row];a<ar[row+1];++a) {
            int shared=ai[a]-1;
            if(shared<0 || shared>=n) std::abort();
            for(int b=br[shared];b<br[shared+1];++b)
                dense.at(bi[b]-1)+=av[a]*bv[b];
        }
        for(int col=0;col<k;++col) if(dense[col]!=0) {
            c.push_back(col); v.push_back(dense[col]);
        }
        r[row+1]=c.size();
    }
    *nz=c.size();
    if(accumulating) {
        if(m!=accumulated_m || k!=accumulated_k) return -1;
        for(int row=0;row<m;++row) for(int j=r[row];j<r[row+1];++j)
            accumulated[(size_t)row*k+c[j]]+=v[j];
        *rp=nullptr; *ci=nullptr; *val=nullptr;
    } else { *rp=r.data(); *ci=c.data(); *val=v.data(); }
    return 0;
}
int cusparse_ksplit_merge_(int ps,int m,int k,int nz,const int *rp,
    const int *ci,const double *v,double thr,int *out_nz,int **out_rp,
    int **out_ci,double **out_v) {
    static std::vector<int> r,c;
    static std::vector<double> values;
    std::vector<double> dense(m*k,0);
    if(accumulating) {
        if(m!=accumulated_m || k!=accumulated_k || nz!=0) return -1;
        dense=accumulated; accumulating=false;
    } else {
    for(int s=0;s<ps;++s) for(int row=0;row<m;++row)
        for(int j=rp[s*(m+1)+row];j<rp[s*(m+1)+row+1];++j) {
            if(j<0 || j>=nz) std::abort();
            dense.at(row*k+ci[j])+=v[j];
        }
    }
    r.assign(m+1,0); c.clear(); values.clear();
    for(int row=0;row<m;++row) {
        for(int col=0;col<k;++col) if(fabs(dense[row*k+col])>thr) {
            c.push_back(col+1); values.push_back(dense[row*k+col]);
        }
        r[row+1]=c.size();
    }
    *out_nz=c.size(); *out_rp=r.data(); *out_ci=c.data(); *out_v=values.data();
    return 0;
}
// Other driver branches must never enter a GPU multiply in this test suite.
int cusparse_gemm_ab_dense_(int m,int n,int k,int na,const int *ar,
    const int *ai,const double *av,int nb,const int *br,const int *bi,
    const double *bv,double threshold,int *nz,int **rp,int **ci,double **v) {
    if(cusparse_solver_begin_(m,k)) return -1;
    int partial,*pr,*pc; double *pv;
    if(cusparse_gemm_ab_partial_(m,n,k,na,ar,ai,av,nb,br,bi,bv,&partial,&pr,&pc,&pv)) return -1;
    int dummy=0; double zero=0;
    return cusparse_ksplit_merge_(1,m,k,0,&dummy,&dummy,&zero,threshold,nz,rp,ci,v);
}
int cusparse_gemm_threshold_() { std::abort(); }
int cusparse_gemm_ab_threshold_() { std::abort(); }
}
