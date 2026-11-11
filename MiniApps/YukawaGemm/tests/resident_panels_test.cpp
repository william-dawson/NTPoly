// Distributed, nonsymmetric CSR-square tests independent of the TRS2 driver.
#include <mpi.h>
#include <cuda_runtime.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cstdio>
extern "C" {
void gpu_ops_init_(int*,int*);
void gpu_ops_alloc_(int*,int*,int*,int**,int**,double**,int*);
void gpu_ops_get_(int*,int**,int**,double**,int*,int*);
int solver_nccl_id_(int8_t*);
int solver_nccl_init_(const int*,const int*,const int*,const int8_t*);
int solver_nccl_members_(const int*,const int*);
int solver_nccl_finish_();
int gpu_resident_square_(int,int,int,int,int,int,int,double,double*);
int gpu_resident_free_();
}
static void check(bool ok) { if(!ok) { std::fprintf(stderr,"Resident CSR test failed\n"); MPI_Abort(MPI_COMM_WORLD,1); } }
int main(int argc,char **argv) {
    MPI_Init(&argc,&argv);
    int rank,size,error,device;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    check(size==1 || size==2 || size==4 || size==8);
    int n=argc>1?std::atoi(argv[1]):6; check(n>0 && n%2==0);
    gpu_ops_init_(&error,&device); check(error==0);
    for(int orientation=0;orientation<(size==2?2:1);++orientation) {
        int pr=size>=4?2:orientation==1?2:1,pc=size/pr;
        check(n%pr==0 && n%pc==0);
        int m=n/pr,k=n/pc,row=rank/pc,col=rank%pc;
        MPI_Comm comms[3]={MPI_COMM_WORLD,MPI_COMM_NULL,MPI_COMM_NULL};
        MPI_Comm_split(MPI_COMM_WORLD,row,col,&comms[1]);
        MPI_Comm_split(MPI_COMM_WORLD,col,row,&comms[2]);
        for(int group=0;group<3;++group) {
            int cr,nr; MPI_Comm_rank(comms[group],&cr); MPI_Comm_size(comms[group],&nr);
            int8_t id[128]={}; const char *shared=std::getenv("YUKAWA_NCCL_SHARED_COMM");
            if(cr==0 && nr>1 && (group==0 || !shared || std::strcmp(shared,"1"))) check(solver_nccl_id_(id)==0);
            MPI_Bcast(id,128,MPI_BYTE,0,comms[group]);
            check(solver_nccl_init_(&group,&nr,&cr,id)==0);
            std::vector<int> members(nr);
            MPI_Allgather(&rank,1,MPI_INT,members.data(),1,MPI_INT,comms[group]);
            check(solver_nccl_members_(&group,members.data())==0);
        }
        for(int example=0;example<3;++example) {
            std::vector<double> matrix(size_t(n)*n,0);
            auto at=[&](int r,int c)->double& { return matrix[size_t(r)*n+c]; };
            for(int r=0;r<n;++r) for(int c=0;c<n;++c)
                if(example!=1 && r!=n-1 && c!=0 && (r*3+c)%4)
                    at(r,c)=double((r*7+c*11)%9-4)/8;
            if(example==1) at(0,0)=1;
            std::vector<int> rp(1,0),ci; std::vector<double> values;
            for(int r=0;r<m;++r) {
                for(int c=0;c<k;++c) if(at(row*m+r,col*k+c)!=0) {
                    ci.push_back(c+1); values.push_back(at(row*m+r,col*k+c));
                }
                rp.push_back(int(ci.size()));
            }
            int slot=1,nz=int(ci.size()),*drp,*dci; double *dv;
            gpu_ops_alloc_(&slot,&m,&nz,&drp,&dci,&dv,&error);check(error==0);
            std::memcpy(drp,rp.data(),sizeof(int)*rp.size());
            if(nz) { std::memcpy(dci,ci.data(),sizeof(int)*nz); std::memcpy(dv,values.data(),sizeof(double)*nz); }
            double phases[4];
            check(gpu_resident_square_(1,2,m,n,k,pr,pc,.125,phases)==0);
            slot=2; int rows;
            gpu_ops_get_(&slot,&drp,&dci,&dv,&rows,&nz);
            check(rows==m && drp[0]==0 && drp[m]==nz);
            int cursor=0;
            for(int r=0;r<m;++r) {
                check(drp[r]==cursor);
                for(int c=0;c<k;++c) {
                    double expected=0;
                    for(int t=0;t<n;++t) expected+=at(row*m+r,t)*at(t,col*k+c);
                    if(std::fabs(expected)>.125) {
                        check(cursor<nz && dci[cursor]==c+1 && std::fabs(dv[cursor]-expected)<1e-14);
                        ++cursor;
                    }
                }
                check(drp[r+1]==cursor);
            }
            check(cursor==nz);
        }
        check(gpu_resident_free_()==0); check(solver_nccl_finish_()==0);
        MPI_Comm_free(&comms[1]); MPI_Comm_free(&comms[2]);
        if(rank==0) std::printf("Resident nonsymmetric/empty/reset tests: PASS grid=%dx%d N=%d\n",pr,pc,n);
    }
    MPI_Finalize();
}
