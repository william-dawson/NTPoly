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
#include <new>
#ifdef YUKAWA_HAVE_NCCL
#include <cuda_runtime.h>
#include <nccl.h>
namespace {
struct Group {
    ncclComm_t comm=nullptr;
    int size=0,rank=0;
    bool world_shared=false;
    std::vector<int> members;
};
struct Resource {
    void* pointer=nullptr;
    void (*destroy)(void*)=nullptr;
    uint64_t (*bytes)(void*)=nullptr;
};
struct WorkspaceOwner {
    Resource resource[3]; // complex, real algebra, CSR panels
    int references=1,device=-1;
    bool leased=false;
};
struct Context {
    Group groups[3];
    cudaStream_t stream=nullptr;
    cudaEvent_t device_ready=nullptr,device_done=nullptr;
    void *input=nullptr,*output=nullptr;
    size_t input_capacity=0,output_capacity=0;
    int device=-1;
    WorkspaceOwner workspace;
    WorkspaceOwner* pool=nullptr;
};
// The mini-app retains its legacy context. Library grids own independent contexts.
Context legacy_context;
Context *current=&legacy_context;
int owned_contexts=0;
std::vector<WorkspaceOwner*> workspace_registry;
#define CUDA(call) do { auto e=(call); if(e!=cudaSuccess) { std::fprintf(stderr,"NCCL staging CUDA: %s\n",cudaGetErrorString(e)); return -1; } } while(0)
#define NCCL(call) do { auto e=(call); if(e!=ncclSuccess) { std::fprintf(stderr,"NCCL: %s\n",ncclGetErrorString(e)); return -1; } } while(0)
WorkspaceOwner* active_workspace(){return current->pool?current->pool:&current->workspace;}
uint64_t workspace_bytes(WorkspaceOwner* owner) {
    uint64_t bytes=0;for(auto& resource:owner->resource)if(resource.pointer)bytes+=resource.bytes(resource.pointer);
    return bytes;
}
int clear_workspace(WorkspaceOwner* owner) {
    CUDA(cudaDeviceSynchronize());
    for(auto& resource:owner->resource)if(resource.pointer){resource.destroy(resource.pointer);resource=Resource{};}
    return 0;
}
int clear_idle_workspace(WorkspaceOwner* owner,int64_t limit,int* cleared) {
    *cleared=0;
    if(limit<0 || owner->leased || (current->groups[0].size>0 && active_workspace()==owner))return 0;
    auto bytes=workspace_bytes(owner);
    if(bytes==0 || (limit>0 && bytes<=uint64_t(limit)))return 0;
    int device;CUDA(cudaGetDevice(&device));CUDA(cudaSetDevice(owner->device));
    int status=clear_workspace(owner);if(status==0)*cleared=1;
    CUDA(cudaSetDevice(device));return status;
}
int reserve(void **p,size_t &capacity,size_t required) {
    required=std::max(required,size_t(8));
    if(capacity>=required) return 0;
    if(*p) { CUDA(cudaFree(*p)); *p=nullptr; capacity=0; }
    CUDA(cudaMalloc(p,required)); capacity=required; return 0;
}
int gather(int which,const void *send,int count,void *receive,const int *counts,
           const int *displacements,bool real) {
    if(which<0 || which>=3 || count<0) return -1;
    Group &g=current->groups[which];
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
    if(reserve(&current->input,current->input_capacity,size_t(count)*width) || reserve(&current->output,current->output_capacity,total*width)) return -1;
    if(count) CUDA(cudaMemcpyAsync(current->input,send,size_t(count)*width,cudaMemcpyHostToDevice,current->stream));
    ncclDataType_t type=real?ncclDouble:ncclInt;
    if(uniform && !g.world_shared) {
        NCCL(ncclAllGather(current->input,current->output,count,type,g.comm,current->stream));
    } else {
        // All peers must make progress together (NCCL grouped P2P).
        NCCL(ncclGroupStart());
        for(int peer=0;peer<g.size;++peer) if(peer!=g.rank) {
            const int target=g.world_shared ? g.members[peer] : peer;
            if(count) NCCL(ncclSend(current->input,count,type,target,g.comm,current->stream));
            if(counts[peer]) NCCL(ncclRecv(static_cast<char*>(current->output)+size_t(displacements[peer])*width,
                                         counts[peer],type,target,g.comm,current->stream));
        }
        NCCL(ncclGroupEnd());
        if(count) CUDA(cudaMemcpyAsync(static_cast<char*>(current->output)+size_t(displacements[g.rank])*width,
                                      current->input,size_t(count)*width,cudaMemcpyDeviceToDevice,current->stream));
    }
    if(total) CUDA(cudaMemcpyAsync(receive,current->output,total*width,cudaMemcpyDeviceToHost,current->stream));
    CUDA(cudaStreamSynchronize(current->stream)); return 0;
}
}
extern "C" {
int solver_nccl_context_create_(void **ptr) {
    if(!ptr || *ptr) return -1;
    auto *context=new(std::nothrow) Context;
    if(!context) return -1;
    auto error=cudaGetDevice(&context->device);
    if(error!=cudaSuccess) { delete context; return -1; }
    context->workspace.device=context->device;
    try {workspace_registry.push_back(&context->workspace);} catch(const std::bad_alloc&){delete context;return -1;}
    *ptr=context;
    ++owned_contexts;
    return 0;
}
int solver_nccl_context_count_() { return owned_contexts; }
int solver_nccl_context_select_(void *ptr) {
    if(!ptr || current!=&legacy_context || current->groups[0].size>0) return -1;
    auto *context=static_cast<Context*>(ptr);
    int device; CUDA(cudaGetDevice(&device));
    if(device!=context->device) return -1;
    current=context;
    return 0;
}
// Grid caches and explicit memory pools share the same typed scratch interface.
int solver_nccl_workspace_get_kind_(int kind,void** pointer) {
    if(!pointer || kind<0 || kind>=3 || current->groups[0].size<1)return -1;
    *pointer=active_workspace()->resource[kind].pointer;return 0;
}
int solver_nccl_workspace_attach_kind_(int kind,void* pointer,void(*destroy)(void*),uint64_t(*bytes)(void*)) {
    if(kind<0 || kind>=3 || !pointer || !destroy || !bytes || current->groups[0].size<1)return -1;
    auto& resource=active_workspace()->resource[kind];if(resource.pointer)return -1;
    resource={pointer,destroy,bytes};return 0;
}
int solver_nccl_workspace_get_(void** pointer){return solver_nccl_workspace_get_kind_(0,pointer);}
int solver_nccl_workspace_attach_(void* p,void(*destroy)(void*),uint64_t(*bytes)(void*)) {
    return solver_nccl_workspace_attach_kind_(0,p,destroy,bytes);
}
int solver_nccl_workspace_clear_kind_(int kind) {
    if(kind<0 || kind>=3)return -1;
    auto& resource=active_workspace()->resource[kind];
    if(resource.pointer){CUDA(cudaDeviceSynchronize());resource.destroy(resource.pointer);resource=Resource{};}
    return 0;
}
int solver_nccl_workspace_clear_(){return clear_workspace(active_workspace());}
int solver_nccl_context_clear_workspace_(void* pointer,const int64_t* limit,int* cleared) {
    *cleared=0;if(!pointer || *limit<0)return 0;
    auto* context=static_cast<Context*>(pointer);if(context==current)return -1;
    return clear_idle_workspace(&context->workspace,*limit,cleared);
}
int solver_nccl_trim_workspaces_limit_(const int64_t* limit,int* cleared) {
    *cleared=0;
    for(auto* owner:workspace_registry){int hit=0;int status=clear_idle_workspace(owner,*limit,&hit);
        if(status)return status;*cleared+=hit;}
    return 0;
}
int solver_nccl_trim_workspaces_(int* cleared) {
    const int64_t limit=0;return solver_nccl_trim_workspaces_limit_(&limit,cleared);
}
int solver_gpu_pool_select_(void** pointer) {
    if(!pointer || current->groups[0].size<1 || current->pool)return -1;
    if(!*pointer){
        auto* owner=new(std::nothrow) WorkspaceOwner;if(!owner)return -1;
        owner->device=current->workspace.device;
        try{workspace_registry.push_back(owner);}catch(const std::bad_alloc&){delete owner;return -1;}
        *pointer=owner;
    }
    auto* owner=static_cast<WorkspaceOwner*>(*pointer);
    int device;CUDA(cudaGetDevice(&device));
    if(owner->device!=device || owner->leased)return -1;
    owner->leased=true;current->pool=owner;return 0;
}
int solver_gpu_pool_retain_(void* pointer) {
    if(!pointer)return 0;auto* owner=static_cast<WorkspaceOwner*>(pointer);
    if(owner->references==INT32_MAX)return -1;++owner->references;return 0;
}
int solver_gpu_pool_release_(void** pointer) {
    if(!pointer || !*pointer)return 0;auto* owner=static_cast<WorkspaceOwner*>(*pointer);
    if(owner->leased)return -1;
    if(owner->references>1){--owner->references;*pointer=nullptr;return 0;}
    int device;CUDA(cudaGetDevice(&device));CUDA(cudaSetDevice(owner->device));
    int status=clear_workspace(owner);CUDA(cudaSetDevice(device));if(status)return status;
    workspace_registry.erase(std::remove(workspace_registry.begin(),workspace_registry.end(),owner),workspace_registry.end());
    delete owner;*pointer=nullptr;return 0;
}
int solver_nccl_workspace_count_(const int* kind) {
    if(*kind<0 || *kind>=3)return -1;int count=0;
    for(auto* owner:workspace_registry)if(owner->resource[*kind].pointer)++count;
    return count;
}
int solver_nccl_context_suspend_() {
    if(current==&legacy_context)return -1;
    CUDA(cudaDeviceSynchronize());
    if(current->pool){current->pool->leased=false;current->pool=nullptr;}
    current=&legacy_context;return 0;
}
int solver_nccl_available_() { return 1; }
int solver_nccl_id_(int8_t *bytes) {
    static_assert(sizeof(ncclUniqueId)==128,"Fortran bootstrap ID size");
    ncclUniqueId id; NCCL(ncclGetUniqueId(&id)); std::memcpy(bytes,&id,128); return 0;
}
int solver_nccl_init_(const int *which,const int *size,const int *rank,const int8_t *bytes) {
    if(*which<0 || *which>=3 || *size<1 || *rank<0 || *rank>=*size) return -1;
    if(current->workspace.device<0){
        CUDA(cudaGetDevice(&current->workspace.device));
        try{workspace_registry.push_back(&current->workspace);}catch(const std::bad_alloc&){return -1;}
    }
    if(!current->stream) CUDA(cudaStreamCreateWithFlags(&current->stream,cudaStreamNonBlocking));
    Group &g=current->groups[*which]; g.size=*size; g.rank=*rank;
    const char *shared=std::getenv("YUKAWA_NCCL_SHARED_COMM");
    g.world_shared=*which>0 && shared && std::strcmp(shared,"1")==0;
    if(*size>1) {
        if(g.world_shared) {
            if(!current->groups[0].comm) return -1;
            g.comm=current->groups[0].comm;
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
    Group &g=current->groups[*which];
    g.members.assign(members,members+g.size);
    for(int rank:g.members) if(rank<0 || rank>=current->groups[0].size) return -1;
    return 0;
}
int solver_nccl_group_info_(int which,int *size,int *rank) {
    if(which<0 || which>=3 || current->groups[which].size<1) return -1;
    *size=current->groups[which].size; *rank=current->groups[which].rank; return 0;
}
// Device-facing transport. Producer/consumer work stays on the default
// compute current->stream; events order it against the nonblocking NCCL current->stream.
// Only counts/displacements live on the CPU. No payload staging or host wait.
int solver_nccl_device_gather_(int which,const void *send,int count,void *receive,
                               const int *counts,const int *displs,int real) {
    if(which<0 || which>=3 || count<0) return -1;
    Group &g=current->groups[which];
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
    if(!current->device_ready) CUDA(cudaEventCreateWithFlags(&current->device_ready,cudaEventDisableTiming));
    if(!current->device_done) CUDA(cudaEventCreateWithFlags(&current->device_done,cudaEventDisableTiming));
    CUDA(cudaEventRecord(current->device_ready,0));
    CUDA(cudaStreamWaitEvent(current->stream,current->device_ready,0));
    const auto type=real?ncclDouble:ncclInt;
    if(uniform && !g.world_shared) {
        NCCL(ncclAllGather(send,receive,count,type,g.comm,current->stream));
    } else {
        NCCL(ncclGroupStart());
        for(int peer=0;peer<g.size;++peer) if(peer!=g.rank) {
            int target=g.world_shared?g.members[peer]:peer;
            if(count) NCCL(ncclSend(send,count,type,target,g.comm,current->stream));
            if(counts[peer]) NCCL(ncclRecv(static_cast<char*>(receive)+size_t(displs[peer])*width,
                                         counts[peer],type,target,g.comm,current->stream));
        }
        NCCL(ncclGroupEnd());
        if(count) CUDA(cudaMemcpyAsync(static_cast<char*>(receive)+size_t(displs[g.rank])*width,
                                      send,size_t(count)*width,cudaMemcpyDeviceToDevice,current->stream));
    }
    CUDA(cudaEventRecord(current->device_done,current->stream));
    CUDA(cudaStreamWaitEvent(0,current->device_done,0));
    return 0;
}
int solver_nccl_gather_i_(const int *which,const int *send,const int *count,int *receive,const int *counts,const int *displacements) {
    return gather(*which,send,*count,receive,counts,displacements,false);
}
int solver_nccl_gather_r_(const int *which,const double *send,const int *count,double *receive,const int *counts,const int *displacements) {
    return gather(*which,send,*count,receive,counts,displacements,true);
}
int solver_nccl_sum_(double *value) {
    Group &g=current->groups[0]; if(g.size==1) return 0; if(g.size<1) return -1;
    if(reserve(&current->input,current->input_capacity,sizeof(double)) || reserve(&current->output,current->output_capacity,sizeof(double))) return -1;
    CUDA(cudaMemcpyAsync(current->input,value,sizeof(double),cudaMemcpyHostToDevice,current->stream));
    NCCL(ncclAllReduce(current->input,current->output,1,ncclDouble,ncclSum,g.comm,current->stream));
    CUDA(cudaMemcpyAsync(value,current->output,sizeof(double),cudaMemcpyDeviceToHost,current->stream));
    CUDA(cudaStreamSynchronize(current->stream)); return 0;
}
int solver_nccl_finish_() {
    if(current->pool)return -1;
    if(clear_workspace(&current->workspace)!=0)return -1;
    if(current->stream) CUDA(cudaStreamSynchronize(current->stream));
    for(auto &g:current->groups) {
        if(g.comm && !g.world_shared) NCCL(ncclCommDestroy(g.comm));
        g.comm=nullptr; g.members.clear(); g.size=0; g.rank=0; g.world_shared=false;
    }
    if(current->input) CUDA(cudaFree(current->input)); if(current->output) CUDA(cudaFree(current->output));
    current->input=current->output=nullptr; current->input_capacity=current->output_capacity=0;
    if(current->device_ready) CUDA(cudaEventDestroy(current->device_ready));
    if(current->device_done) CUDA(cudaEventDestroy(current->device_done));
    current->device_ready=current->device_done=nullptr;
    if(current->stream) CUDA(cudaStreamDestroy(current->stream)); current->stream=nullptr; return 0;
}
int solver_nccl_context_free_(void **ptr) {
    if(!ptr || !*ptr) return 0;
    auto *context=static_cast<Context*>(*ptr);
    if(current==context) return -1; // A leased grid cannot be destroyed.
    int old_device; CUDA(cudaGetDevice(&old_device));
    CUDA(cudaSetDevice(context->device));
    Context *saved=current; current=context;
    int status=solver_nccl_finish_(); current=saved;
    CUDA(cudaSetDevice(old_device));
    if(status) return status;
    workspace_registry.erase(std::remove(workspace_registry.begin(),workspace_registry.end(),&context->workspace),workspace_registry.end());
    delete context; *ptr=nullptr; --owned_contexts;
    return 0;
}
}
#else
extern "C" {
int solver_nccl_context_count_() { return 0; }
int solver_nccl_workspace_get_kind_(int,void**){return -1;}
int solver_nccl_workspace_attach_kind_(int,void*,void(*)(void*),uint64_t(*)(void*)){return -1;}
int solver_nccl_workspace_clear_kind_(int){return 0;}
int solver_nccl_trim_workspaces_limit_(const int64_t*,int*){return -1;}
int solver_gpu_pool_select_(void**){return -1;}
int solver_gpu_pool_retain_(void*){return -1;}
int solver_gpu_pool_release_(void**){return -1;}
int solver_nccl_workspace_count_(const int*){return 0;}

int solver_nccl_workspace_get_(void**) { return -1; }
int solver_nccl_workspace_attach_(void*,void(*)(void*),uint64_t(*)(void*)) { return -1; }
int solver_nccl_workspace_clear_() { return 0; }
int solver_nccl_trim_workspaces_(int*) { return -1; }
int solver_nccl_context_clear_workspace_(void*,const int64_t*,int*) { return -1; }
int solver_nccl_context_create_(void**) { return -1; }
int solver_nccl_context_select_(void*) { return -1; }
int solver_nccl_context_suspend_() { return -1; }
int solver_nccl_context_free_(void**) { return -1; }
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
