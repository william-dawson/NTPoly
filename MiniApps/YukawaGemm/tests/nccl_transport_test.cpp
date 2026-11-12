#include <mpi.h>
#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cstring>
extern "C" {
int solver_nccl_available_();
int solver_nccl_workspace_get_(void**);
int solver_nccl_workspace_attach_(void*,void(*)(void*),uint64_t(*)(void*));
int solver_nccl_id_(int8_t*);
int solver_nccl_init_(const int*,const int*,const int*,const int8_t*);
int solver_nccl_members_(const int*,const int*);
int solver_nccl_gather_i_(const int*,const int*,const int*,int*,const int*,const int*);
int solver_nccl_gather_r_(const int*,const double*,const int*,double*,const int*,const int*);
int solver_nccl_sum_(double*);
int solver_nccl_finish_();
}
static int destroyed=0;
static void release_workspace(void* p){delete static_cast<int*>(p);++destroyed;}
static uint64_t workspace_bytes(void*){return sizeof(int);}
static void check(bool ok) { if(!ok) MPI_Abort(MPI_COMM_WORLD,1); }
int main(int argc,char **argv) {
    MPI_Init(&argc,&argv);
    int rank,size; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    check(size==1 || size==2 || size==4 || size==8);
    MPI_Comm local;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
    int local_rank; MPI_Comm_rank(local,&local_rank);
    int devices=0; check(cudaGetDeviceCount(&devices)==cudaSuccess && devices>0);
    check(cudaSetDevice(local_rank%devices)==cudaSuccess && solver_nccl_available_());
    MPI_Comm_free(&local);
    int columns=size>=4?size/2:size;
    MPI_Comm comms[3]={MPI_COMM_WORLD,MPI_COMM_NULL,MPI_COMM_NULL};
    MPI_Comm_split(MPI_COMM_WORLD,rank/columns,rank,&comms[1]);
    MPI_Comm_split(MPI_COMM_WORLD,rank%columns,rank,&comms[2]);
    for(int group=0;group<3;++group) {
        int cr,nr; MPI_Comm_rank(comms[group],&cr); MPI_Comm_size(comms[group],&nr);
        int8_t id[128]={};
        const char *shared=std::getenv("YUKAWA_NCCL_SHARED_COMM");
        if(cr==0 && nr>1 && (group==0 || !shared || std::strcmp(shared,"1"))) check(solver_nccl_id_(id)==0);
        MPI_Bcast(id,128,MPI_BYTE,0,comms[group]);
        check(solver_nccl_init_(&group,&nr,&cr,id)==0);
        std::vector<int> members(nr);
        MPI_Allgather(&rank,1,MPI_INT,members.data(),1,MPI_INT,comms[group]);
        check(solver_nccl_members_(&group,members.data())==0);
    }
    void* stored=nullptr;
    check(solver_nccl_workspace_get_(&stored)==0 && !stored);
    auto* resource=new int(7);
    check(solver_nccl_workspace_attach_(resource,release_workspace,workspace_bytes)==0);
    check(solver_nccl_workspace_attach_(resource,release_workspace,workspace_bytes)!=0);
    check(solver_nccl_workspace_get_(&stored)==0 && stored==resource && *static_cast<int*>(stored)==7);
    for(int group=0;group<3;++group) {
        int cr,nr; MPI_Comm_rank(comms[group],&cr); MPI_Comm_size(comms[group],&nr);
        std::vector<int> counts(nr,2),displs(nr),received(2*nr);
        for(int p=0;p<nr;++p) displs[p]=2*p;
        int count=2,send[2]={cr,100+cr};
        check(solver_nccl_gather_i_(&group,send,&count,received.data(),counts.data(),displs.data())==0);
        for(int p=0;p<nr;++p) check(received[2*p]==p && received[2*p+1]==100+p);
        // Variable counts, including an empty contributor and all-empty gathers.
        int total=0;
        for(int p=0;p<nr;++p) { counts[p]=p==0?0:p+1; displs[p]=total; total+=counts[p]; }
        count=counts[cr];
        std::vector<double> values(count,cr+0.25),result(total);
        check(solver_nccl_gather_r_(&group,values.data(),&count,result.data(),counts.data(),displs.data())==0);
        for(int p=0;p<nr;++p) for(int j=0;j<counts[p];++j) check(result[displs[p]+j]==p+0.25);
        for(int p=0;p<nr;++p) counts[p]=displs[p]=0;
        count=0;
        check(solver_nccl_gather_r_(&group,values.data(),&count,result.data(),counts.data(),displs.data())==0);
    }
    double sum=rank+1;
    check(solver_nccl_sum_(&sum)==0 && sum==size*(size+1)/2.0);
    check(solver_nccl_finish_()==0 && destroyed==1);
    check(solver_nccl_workspace_get_(&stored)!=0);
    if(rank==0) std::printf("NCCL transport fixed/variable/empty/scalar checks: PASS np=%d\n",size);
    MPI_Comm_free(&comms[1]); MPI_Comm_free(&comms[2]);
    MPI_Finalize();
}
