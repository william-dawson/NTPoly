// Blocking host-facing adapters around device-buffer NCCL communication.
// The solver still constructs CSR panels on the host. Each adapter stages
// the whole payload in one transfer in each direction, never via MPI's
// small CUDA-aware copies into managed receive arrays.
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>
#ifdef YUKAWA_HAVE_NCCL
#include <cuda_runtime.h>
#include <nccl.h>
namespace {
struct Group {
    ncclComm_t comm=nullptr;
    int size=0,rank=0;
    bool world_shared=false;
    std::vector<int> members;
} groups[3];
cudaStream_t stream=nullptr;
cudaEvent_t device_ready=nullptr,device_done=nullptr;
void *input=nullptr,*output=nullptr;
size_t input_capacity=0,output_capacity=0;
#define CUDA(call) do { auto e=(call); if(e!=cudaSuccess) { std::fprintf(stderr,"NCCL staging CUDA: %s\n",cudaGetErrorString(e)); return -1; } } while(0)
#define NCCL(call) do { auto e=(call); if(e!=ncclSuccess) { std::fprintf(stderr,"NCCL: %s\n",ncclGetErrorString(e)); return -1; } } while(0)
int reserve(void **p,size_t &capacity,size_t required) {
    required=std::max(required,size_t(8));
    if(capacity>=required) return 0;
    if(*p) { CUDA(cudaFree(*p)); *p=nullptr; capacity=0; }
    CUDA(cudaMalloc(p,required)); capacity=required; return 0;
}
int gather(int which,const void *send,int count,void *receive,const int *counts,
           const int *displacements,bool real) {
    if(which<0 || which>=3 || count<0) return -1;
    Group &g=groups[which];
    if(g.size<1 || counts[g.rank]!=count) return -1;
    const size_t width=real?sizeof(double):sizeof(int);
    size_t total=0; bool uniform=true;
    for(int peer=0;peer<g.size;++peer) {
        if(counts[peer]<0 || displacements[peer]<0) return -1;
        total=std::max(total,size_t(displacements[peer])+size_t(counts[peer]));
        uniform=uniform && counts[peer]==count && size_t(displacements[peer])==size_t(peer)*count;
    }
    if(total==0) return 0;
    if(g.size==1) {
        if(count) std::memmove(static_cast<char*>(receive)+displacements[0]*width,send,count*width);
        return 0;
    }
    if(g.world_shared && g.members.size()!=size_t(g.size)) return -1;
    if(reserve(&input,input_capacity,size_t(count)*width) || reserve(&output,output_capacity,total*width)) return -1;
    if(count) CUDA(cudaMemcpyAsync(input,send,size_t(count)*width,cudaMemcpyHostToDevice,stream));
    ncclDataType_t type=real?ncclDouble:ncclInt;
    if(uniform && !g.world_shared) {
        NCCL(ncclAllGather(input,output,count,type,g.comm,stream));
    } else {
        // All peers must make progress together (NCCL grouped P2P).
        NCCL(ncclGroupStart());
        for(int peer=0;peer<g.size;++peer) if(peer!=g.rank) {
            const int target=g.world_shared ? g.members[peer] : peer;
            if(count) NCCL(ncclSend(input,count,type,target,g.comm,stream));
            if(counts[peer]) NCCL(ncclRecv(static_cast<char*>(output)+size_t(displacements[peer])*width,
                                         counts[peer],type,target,g.comm,stream));
        }
        NCCL(ncclGroupEnd());
        if(count) CUDA(cudaMemcpyAsync(static_cast<char*>(output)+size_t(displacements[g.rank])*width,
                                      input,size_t(count)*width,cudaMemcpyDeviceToDevice,stream));
    }
    if(total) CUDA(cudaMemcpyAsync(receive,output,total*width,cudaMemcpyDeviceToHost,stream));
    CUDA(cudaStreamSynchronize(stream)); return 0;
}
}
extern "C" {
int solver_nccl_available_() { return 1; }
int solver_nccl_id_(int8_t *bytes) {
    static_assert(sizeof(ncclUniqueId)==128,"Fortran bootstrap ID size");
    ncclUniqueId id; NCCL(ncclGetUniqueId(&id)); std::memcpy(bytes,&id,128); return 0;
}
int solver_nccl_init_(const int *which,const int *size,const int *rank,const int8_t *bytes) {
    if(*which<0 || *which>=3 || *size<1 || *rank<0 || *rank>=*size) return -1;
    if(!stream) CUDA(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    Group &g=groups[*which]; g.size=*size; g.rank=*rank;
    const char *shared=std::getenv("YUKAWA_NCCL_SHARED_COMM");
    g.world_shared=*which>0 && shared && std::strcmp(shared,"1")==0;
    if(*size>1) {
        if(g.world_shared) {
            if(!groups[0].comm) return -1;
            g.comm=groups[0].comm;
        } else {
            ncclUniqueId id; std::memcpy(&id,bytes,128);
            NCCL(ncclCommInitRank(&g.comm,*size,id,*rank));
        }
    }
    if(*which==0 && *rank==0) { int version; NCCL(ncclGetVersion(&version)); std::printf("NCCL version: %d\n",version); }
    return 0;
}
int solver_nccl_members_(const int *which,const int *members) {
    if(*which<0 || *which>=3) return -1;
    Group &g=groups[*which];
    g.members.assign(members,members+g.size);
    for(int rank:g.members) if(rank<0 || rank>=groups[0].size) return -1;
    return 0;
}
int solver_nccl_group_info_(int which,int *size,int *rank) {
    if(which<0 || which>=3 || groups[which].size<1) return -1;
    *size=groups[which].size; *rank=groups[which].rank; return 0;
}
// Device-facing transport. Producer/consumer work stays on the default
// compute stream; events order it against the nonblocking NCCL stream.
// Only counts/displacements live on the CPU. No payload staging or host wait.
int solver_nccl_device_gather_(int which,const void *send,int count,void *receive,
                               const int *counts,const int *displs,int real) {
    if(which<0 || which>=3 || count<0) return -1;
    Group &g=groups[which];
    if(g.size<1 || counts[g.rank]!=count) return -1;
    const size_t width=real?sizeof(double):sizeof(int);
    size_t total=0; bool uniform=true;
    for(int p=0;p<g.size;++p) {
        if(counts[p]<0 || displs[p]<0) return -1;
        total=std::max(total,size_t(displs[p])+size_t(counts[p]));
        uniform=uniform && counts[p]==count && size_t(displs[p])==size_t(p)*count;
    }
    if(!total) return 0;
    if(g.size==1) {
        if(count) CUDA(cudaMemcpyAsync(static_cast<char*>(receive)+size_t(displs[0])*width,
                                      send,size_t(count)*width,cudaMemcpyDeviceToDevice));
        return 0;
    }
    if(g.world_shared && g.members.size()!=size_t(g.size)) return -1;
    if(!device_ready) CUDA(cudaEventCreateWithFlags(&device_ready,cudaEventDisableTiming));
    if(!device_done) CUDA(cudaEventCreateWithFlags(&device_done,cudaEventDisableTiming));
    CUDA(cudaEventRecord(device_ready,0));
    CUDA(cudaStreamWaitEvent(stream,device_ready,0));
    const auto type=real?ncclDouble:ncclInt;
    if(uniform && !g.world_shared) {
        NCCL(ncclAllGather(send,receive,count,type,g.comm,stream));
    } else {
        NCCL(ncclGroupStart());
        for(int peer=0;peer<g.size;++peer) if(peer!=g.rank) {
            int target=g.world_shared?g.members[peer]:peer;
            if(count) NCCL(ncclSend(send,count,type,target,g.comm,stream));
            if(counts[peer]) NCCL(ncclRecv(static_cast<char*>(receive)+size_t(displs[peer])*width,
                                         counts[peer],type,target,g.comm,stream));
        }
        NCCL(ncclGroupEnd());
        if(count) CUDA(cudaMemcpyAsync(static_cast<char*>(receive)+size_t(displs[g.rank])*width,
                                      send,size_t(count)*width,cudaMemcpyDeviceToDevice,stream));
    }
    CUDA(cudaEventRecord(device_done,stream));
    CUDA(cudaStreamWaitEvent(0,device_done,0));
    return 0;
}
int solver_nccl_gather_i_(const int *which,const int *send,const int *count,int *receive,const int *counts,const int *displacements) {
    return gather(*which,send,*count,receive,counts,displacements,false);
}
int solver_nccl_gather_r_(const int *which,const double *send,const int *count,double *receive,const int *counts,const int *displacements) {
    return gather(*which,send,*count,receive,counts,displacements,true);
}
int solver_nccl_sum_(double *value) {
    Group &g=groups[0]; if(g.size==1) return 0; if(g.size<1) return -1;
    if(reserve(&input,input_capacity,sizeof(double)) || reserve(&output,output_capacity,sizeof(double))) return -1;
    CUDA(cudaMemcpyAsync(input,value,sizeof(double),cudaMemcpyHostToDevice,stream));
    NCCL(ncclAllReduce(input,output,1,ncclDouble,ncclSum,g.comm,stream));
    CUDA(cudaMemcpyAsync(value,output,sizeof(double),cudaMemcpyDeviceToHost,stream));
    CUDA(cudaStreamSynchronize(stream)); return 0;
}
int solver_nccl_finish_() {
    if(stream) CUDA(cudaStreamSynchronize(stream));
    for(auto &g:groups) {
        if(g.comm && !g.world_shared) NCCL(ncclCommDestroy(g.comm));
        g.comm=nullptr; g.members.clear();
    }
    if(input) CUDA(cudaFree(input)); if(output) CUDA(cudaFree(output));
    input=output=nullptr; input_capacity=output_capacity=0;
    if(device_ready) CUDA(cudaEventDestroy(device_ready));
    if(device_done) CUDA(cudaEventDestroy(device_done));
    device_ready=device_done=nullptr;
    if(stream) CUDA(cudaStreamDestroy(stream)); stream=nullptr; return 0;
}
}
#else
extern "C" {
int solver_nccl_available_() { return 0; }
int solver_nccl_id_(int8_t*) { return -1; }
int solver_nccl_init_(const int*,const int*,const int*,const int8_t*) { return -1; }
int solver_nccl_members_(const int*,const int*) { return -1; }
int solver_nccl_group_info_(int,int*,int*) { return -1; }
int solver_nccl_device_gather_(int,const void*,int,void*,const int*,const int*,int) { return -1; }
int solver_nccl_gather_i_(const int*,const int*,const int*,int*,const int*,const int*) { return -1; }
int solver_nccl_gather_r_(const int*,const double*,const int*,double*,const int*,const int*) { return -1; }
int solver_nccl_sum_(double*) { return -1; }
int solver_nccl_finish_() { return 0; }
}
#endif
