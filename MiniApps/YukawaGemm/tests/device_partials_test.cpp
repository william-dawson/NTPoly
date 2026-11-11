// Real CUDA/cuSPARSE regression: sum before threshold, cancellation, empties,
// repeated accumulator reset, and the no-host-partial pointer contract.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
extern "C" {
int cusparse_solver_begin_(int,int);
int cusparse_gemm_ab_dense_(int,int,int,int,const int*,const int*,const double*,
    int,const int*,const int*,const double*,double,int*,int**,int**,double**);
int cusparse_gemm_ab_partial_(int,int,int,int,const int*,const int*,const double*,
    int,const int*,const int*,const double*,int*,int**,int**,double**);
int cusparse_ksplit_merge_(int,int,int,int,const int*,const int*,const double*,
    double,int*,int**,int**,double**);
}
static void check(bool ok) { if(!ok) { std::fprintf(stderr,"device partial test failed\n"); std::exit(1); } }
int main() {
    const int ar[]={0,1,2}, ai[]={1,1}, br[]={0,2}, bi[]={1,2};
    const double av[]={1,1};
    for(int repeat=0;repeat<3;++repeat) {
        check(cusparse_solver_begin_(2,2)==0);
        check(cusparse_solver_begin_(2,2)!=0); // Accidental reset loses work.
        for(int chunk=0;chunk<2;++chunk) {
            const double bv[]={.06,chunk ? -.2 : .2};
            int nz=-1,*rp=nullptr,*ci=nullptr; double *v=nullptr;
            check(cusparse_gemm_ab_partial_(2,1,2,repeat==1 ? 0 : 2,ar,ai,av,
                  2,br,bi,bv,&nz,&rp,&ci,&v)==0);
            check(nz==(repeat==1 ? 0 : 4) && !rp && !ci && !v);
        }
        int nz=-1,*rp=nullptr,*ci=nullptr; double *v=nullptr;
        const int dummy=0; const double zero=0;
        check(cusparse_ksplit_merge_(2,2,2,0,&dummy,&dummy,&zero,.1,&nz,&rp,&ci,&v)==0);
        check(rp && rp[0]==0);
        if(repeat==1) check(nz==0 && rp[1]==0 && rp[2]==0);
        else {
            check(nz==2 && rp[1]==1 && rp[2]==2);
            for(int j=0;j<2;++j) check(ci[j]==1 && std::fabs(v[j]-.12)<1e-14);
        }
    }
    // Cross several warp tiles and test a non-multiple-of-32 column count.
    const int wide_br[]={0,67};
    std::vector<int> wide_bi(67);
    std::vector<double> wide_bv(67);
    for(int c=0;c<67;++c) { wide_bi[c]=c+1; wide_bv[c]=c%3==0 ? .2 : .03; }
    check(cusparse_solver_begin_(2,67)==0);
    for(int chunk=0;chunk<2;++chunk) {
        int nz,*rp,*ci; double *v;
        check(cusparse_gemm_ab_partial_(2,1,67,2,ar,ai,av,67,wide_br,
            wide_bi.data(),wide_bv.data(),&nz,&rp,&ci,&v)==0);
        check(nz==134 && !rp && !ci && !v);
    }
    int nz,*rp,*ci; double *v; const int dummy=0; const double zero=0;
    check(cusparse_ksplit_merge_(2,2,67,0,&dummy,&dummy,&zero,.1,&nz,&rp,&ci,&v)==0);
    check(nz==46 && rp[0]==0 && rp[1]==23 && rp[2]==46);
    for(int row=0;row<2;++row) for(int j=0;j<23;++j)
        check(ci[row*23+j]==j*3+1 && std::fabs(v[row*23+j]-.4)<1e-14);
    const int dar[]={0,2,4}, dai[]={1,2,1,2}, dbr[]={0,2,4}, dbi[]={1,2,1,2};
    const double dav[]={1,1,1,1}, dbv[]={.06,.2,.06,-.2};
    for(int repeat=0;repeat<3;++repeat) {
        check(cusparse_gemm_ab_dense_(2,2,2,repeat==1 ? 0 : 4,dar,dai,dav,
            4,dbr,dbi,dbv,.1,&nz,&rp,&ci,&v)==0);
        if(repeat==1) check(nz==0 && rp[2]==0);
        else check(nz==2 && rp[1]==1 && rp[2]==2 && ci[0]==1 && ci[1]==1 &&
            std::fabs(v[0]-.12)<1e-14 && std::fabs(v[1]-.12)<1e-14);
    }
    std::puts("Device partial accumulation: threshold/cancellation/empty/reset PASS");
}
