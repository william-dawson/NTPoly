// Native complex FP64 shared products. CPU coordinates are an API boundary;
// row/column panel payloads travel GPU-to-GPU through the leased NCCL context.
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <chrono>
#include <cuComplex.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <exception>
#include <cmath>
#include <cstdio>
#include <vector>
#include <memory>
#include <stdexcept>
extern "C" {
int solver_nccl_workspace_get_(void**);
int solver_nccl_workspace_attach_(void*,void(*)(void*),uint64_t(*)(void*));
int solver_nccl_group_info_(int,int*,int*);
int solver_nccl_gather_i_(const int*,const int*,const int*,int*,const int*,const int*);
int solver_nccl_device_gather_(int,const void*,int,void*,const int*,const int*,int);
}
namespace {
#define CUDA(call) do { auto e=(call);if(e!=cudaSuccess){std::fprintf(stderr,"Complex product: %s: %s\n",#call,cudaGetErrorString(e));return -1;} }while(0)
#define CHECK(call) do {if((call)!=0)return -1;}while(0)
struct Coordinate { int row,column;cuDoubleComplex value; };
static_assert(sizeof(Coordinate)%sizeof(int)==0 && sizeof(Coordinate)<=32,"coordinate wire/planner size");
constexpr int words=sizeof(Coordinate)/sizeof(int);
uint64_t allocation_count=0;
int workspace_count=0;
struct Buffer {
 void* p=nullptr;size_t capacity=0;
 Buffer()=default;
 Buffer(const Buffer&)=delete;Buffer& operator=(const Buffer&)=delete;
 ~Buffer(){if(p)cudaFree(p);}
 int allocate(size_t bytes) {
  bytes=std::max<size_t>(bytes,1);
  if(bytes<=capacity)return 0;
  // Free before growth, rather than transiently retaining two full buffers.
  if(p){CUDA(cudaFree(p));p=nullptr;capacity=0;}
  CUDA(cudaMalloc(&p,bytes));capacity=bytes;++allocation_count;return 0;
 }
 template<class T>T* as(){return static_cast<T*>(p);}
};
struct ScratchPool {
 struct Block {std::unique_ptr<Buffer> storage;bool used=false;};
 std::vector<Block> blocks;
 char* acquire(size_t bytes) {
  Block* selected=nullptr;
  for(auto& block:blocks)if(!block.used&&block.storage->capacity>=bytes)
   if(!selected||block.storage->capacity<selected->storage->capacity)selected=&block;
  if(!selected) {
   for(auto& block:blocks)if(!block.used){selected=&block;break;}
   if(!selected){blocks.push_back({std::make_unique<Buffer>(),false});selected=&blocks.back();}
  }
  if(selected->storage->allocate(bytes)!=0)throw std::bad_alloc();
  selected->used=true;return selected->storage->as<char>();
 }
 void release(char* ptr) {
  for(auto& block:blocks)if(block.storage->p==ptr&&block.used){block.used=false;return;}
  throw std::runtime_error("invalid complex scratch release");
 }
};
struct ScratchAllocator {
 using value_type=char;
 ScratchPool* pool;
 char* allocate(std::ptrdiff_t bytes){if(bytes<0)throw std::bad_alloc();return pool->acquire(size_t(bytes));}
 void deallocate(char* ptr,size_t){pool->release(ptr);}
};
struct Workspace {
 Buffer local_a,local_b,panel_a,panel_b,local_c,arp,brp,dense,counts,packed;
 ScratchPool scratch;
 std::vector<Coordinate> host;
 Workspace(){++workspace_count;}
 ~Workspace(){--workspace_count;}
 uint64_t bytes() const {
  uint64_t total=0;
  for(auto* buffer:{&local_a,&local_b,&panel_a,&panel_b,&local_c,&arp,&brp,&dense,&counts,&packed})total+=buffer->capacity;
  for(const auto& block:scratch.blocks)total+=block.storage->capacity;
  return total;
 }
};
void destroy_workspace(void* ptr){delete static_cast<Workspace*>(ptr);}
uint64_t workspace_bytes(void* ptr){return static_cast<Workspace*>(ptr)->bytes();}
Workspace& get_workspace() {
 void* ptr=nullptr;
 if(solver_nccl_workspace_get_(&ptr)!=0)throw std::runtime_error("complex workspace requires an active context");
 if(!ptr) {
  auto owner=std::make_unique<Workspace>();
  if(solver_nccl_workspace_attach_(owner.get(),destroy_workspace,workspace_bytes)!=0)throw std::runtime_error("complex workspace ownership");
  ptr=owner.release();
 }
 return *static_cast<Workspace*>(ptr);
}
struct Less {
 __host__ __device__ bool operator()(Coordinate a,Coordinate b) const {
  return a.row<b.row||(a.row==b.row&&a.column<b.column);
 }
};
int upload(int count,const int* rows,const int* cols,const cuDoubleComplex* values,Buffer& output,std::vector<Coordinate>& host) {
 host.resize(count);
 for(int i=0;i<count;i++)host[i]={rows[i]-1,cols[i]-1,values[i]};
 CHECK(output.allocate(size_t(count)*sizeof(Coordinate)));
 CUDA(cudaMemcpy(output.p,host.data(),host.size()*sizeof(Coordinate),cudaMemcpyHostToDevice));
 return 0;
}
int exchange(int which,int count,Buffer& input,Buffer& output,int& total) {
 int size,rank;CHECK(solver_nccl_group_info_(which,&size,&rank));
 std::vector<int> counts(size),displs(size),ones(size,1),offsets(size);
 for(int i=0;i<size;i++)offsets[i]=i;
 int one=1;CHECK(solver_nccl_gather_i_(&which,&count,&one,counts.data(),ones.data(),offsets.data()));
 int64_t sum=0;
 for(int i=0;i<size;i++) {
  if(counts[i]<0||sum+counts[i]>INT_MAX/words)return -1;
  displs[i]=int(sum)*words;sum+=counts[i];counts[i]*=words;
 }
 total=int(sum);CHECK(output.allocate(size_t(total)*sizeof(Coordinate)));
 CHECK(solver_nccl_device_gather_(which,input.p,count*words,output.p,counts.data(),displs.data(),0));
 return 0;
}
__global__ void row_counts(const Coordinate* input,int count,int offset,int* rp) {
 for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=blockDim.x*gridDim.x)
  atomicAdd(rp+input[i].row-offset+1,1);
}
__global__ void multiply(int columns,const int* arp,const Coordinate* a,const int* brp,
 const Coordinate* b,int column_offset,cuDoubleComplex* dense) {
 int row=blockIdx.x,lane=threadIdx.x&31,warp=threadIdx.x>>5;
 for(int i=arp[row]+warp;i<arp[row+1];i+=blockDim.x/32) {
  Coordinate x=a[i];
  for(int j=brp[x.column]+lane;j<brp[x.column+1];j+=32) {
   Coordinate y=b[j];cuDoubleComplex value=cuCmul(x.value,y.value);
   auto* dst=dense+size_t(row)*columns+y.column-column_offset;
   atomicAdd(&dst->x,value.x);atomicAdd(&dst->y,value.y);
  }
 }
}
__global__ void prune(cuDoubleComplex* dense,size_t cells,double alpha,double threshold) {
 for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<cells;i+=size_t(blockDim.x)*gridDim.x)
  {
   dense[i].x*=alpha;dense[i].y*=alpha;
   if(hypot(dense[i].x,dense[i].y)<=threshold)dense[i]=make_cuDoubleComplex(0,0);
  }
}
__global__ void increment(cuDoubleComplex* dense,int columns,int row_offset,int column_offset,
 const Coordinate* c,int count,double beta) {
 for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
  Coordinate x=c[i];auto* dst=dense+size_t(x.row-row_offset)*columns+x.column-column_offset;
  dst->x+=beta*x.value.x;dst->y+=beta*x.value.y;
 }
}
__global__ void count_output(const cuDoubleComplex* dense,int rows,int columns,int* counts) {
 int row=blockIdx.x;int local=0;
 for(int j=threadIdx.x;j<columns;j+=blockDim.x) {
  auto v=dense[size_t(row)*columns+j];if((v.x!=0||v.y!=0))++local;
 }
 __shared__ int count[256];count[threadIdx.x]=local;__syncthreads();
 for(int stride=128;stride;stride/=2){if(threadIdx.x<stride)count[threadIdx.x]+=count[threadIdx.x+stride];__syncthreads();}
 if(threadIdx.x==0)counts[row+1]=count[0];
}
__global__ void pack_output(const cuDoubleComplex* dense,int columns,int row_offset,int column_offset,
 const int* rp,Coordinate* output) {
 int row=blockIdx.x,lane=threadIdx.x;int out=rp[row];
 for(int base=0;base<columns;base+=32) {
  int col=base+lane;auto v=col<columns?dense[size_t(row)*columns+col]:make_cuDoubleComplex(0,0);
  bool keep=col<columns&&(v.x!=0||v.y!=0);
  unsigned mask=__ballot_sync(0xffffffffu,keep);
  if(keep)output[out+__popc(mask&((1u<<lane)-1))]={row+row_offset+1,col+column_offset+1,v};
  out+=__popc(mask);
 }
}
double elapsed(std::chrono::steady_clock::time_point start) {
 return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
int compute_product(int* rows,int* columns,int* n,int* row_offset,int* column_offset,
 int* na,Buffer& local_a,int* nb,Buffer& local_b,int* nc,Buffer& local_c,
 double* alpha,double* beta,double* threshold,int* count,Workspace& workspace,double* phases) {
 auto& panel_a=workspace.panel_a;auto& panel_b=workspace.panel_b;
 auto& arp=workspace.arp;auto& brp=workspace.brp;auto& dense=workspace.dense;
 auto& counts=workspace.counts;auto& packed=workspace.packed;
 ScratchAllocator allocator{&workspace.scratch};
 auto policy=thrust::cuda::par(allocator);
  auto start=std::chrono::steady_clock::now();
  int count_a,count_b;
  CHECK(exchange(1,*na,local_a,panel_a,count_a));CHECK(exchange(2,*nb,local_b,panel_b,count_b));
  if(phases){phases[0]=elapsed(start);start=std::chrono::steady_clock::now();}
  // Sort panels and build CSR offsets entirely on the device.
  thrust::sort(policy,panel_a.as<Coordinate>(),panel_a.as<Coordinate>()+count_a,Less{});
  thrust::sort(policy,panel_b.as<Coordinate>(),panel_b.as<Coordinate>()+count_b,Less{});
  CHECK(arp.allocate(size_t(*rows+1)*sizeof(int)));CHECK(brp.allocate(size_t(*n+1)*sizeof(int)));
  CUDA(cudaMemset(arp.p,0,size_t(*rows+1)*sizeof(int)));CUDA(cudaMemset(brp.p,0,size_t(*n+1)*sizeof(int)));
  row_counts<<<256,256>>>(panel_a.as<Coordinate>(),count_a,*row_offset,arp.as<int>());
  row_counts<<<256,256>>>(panel_b.as<Coordinate>(),count_b,0,brp.as<int>());CUDA(cudaGetLastError());
  thrust::inclusive_scan(policy,arp.as<int>(),arp.as<int>()+*rows+1,arp.as<int>());
  thrust::inclusive_scan(policy,brp.as<int>(),brp.as<int>()+*n+1,brp.as<int>());
  if(phases){CUDA(cudaDeviceSynchronize());
   phases[1]=elapsed(start);start=std::chrono::steady_clock::now();}
  size_t cells=size_t(*rows)*size_t(*columns);CHECK(dense.allocate(cells*sizeof(cuDoubleComplex)));
  CUDA(cudaMemset(dense.p,0,cells*sizeof(cuDoubleComplex)));
  multiply<<<*rows,128>>>(*columns,arp.as<int>(),panel_a.as<Coordinate>(),brp.as<int>(),panel_b.as<Coordinate>(),
   *column_offset,dense.as<cuDoubleComplex>());CUDA(cudaGetLastError());
  prune<<<256,256>>>(dense.as<cuDoubleComplex>(),cells,*alpha,*threshold);CUDA(cudaGetLastError());
  if(*nc) {
   increment<<<256,256>>>(dense.as<cuDoubleComplex>(),*columns,*row_offset,*column_offset,local_c.as<Coordinate>(),*nc,*beta);
   CUDA(cudaGetLastError());
  }
  if(phases){CUDA(cudaDeviceSynchronize());
   phases[2]=elapsed(start);start=std::chrono::steady_clock::now();}
  CHECK(counts.allocate(size_t(*rows+1)*sizeof(int)));CUDA(cudaMemset(counts.p,0,sizeof(int)));
  count_output<<<*rows,256>>>(dense.as<cuDoubleComplex>(),*rows,*columns,counts.as<int>());CUDA(cudaGetLastError());
  thrust::inclusive_scan(policy,counts.as<int>(),counts.as<int>()+*rows+1,counts.as<int>());
  CUDA(cudaMemcpy(count,counts.as<int>()+*rows,sizeof(int),cudaMemcpyDeviceToHost));
  CHECK(packed.allocate(size_t(*count)*sizeof(Coordinate)));
  pack_output<<<*rows,32>>>(dense.as<cuDoubleComplex>(),*columns,*row_offset,*column_offset,counts.as<int>(),packed.as<Coordinate>());
  CUDA(cudaGetLastError());
  if(phases){CUDA(cudaDeviceSynchronize());phases[3]=elapsed(start);}
  return 0;
}
struct Result { std::vector<Coordinate> entries; };
}
extern "C" int ntpoly_complex_product_(int* rows,int* columns,int* n,int* row_offset,int* column_offset,
 int* na,const int* ar,const int* ac,const cuDoubleComplex* av,
 int* nb,const int* br,const int* bc,const cuDoubleComplex* bv,
 int* nc,const int* cr,const int* cc,const cuDoubleComplex* cv,
 double* alpha,double* beta,double* threshold,void** output,int* count) {
 *output=nullptr;*count=0;
 try {
  auto& workspace=get_workspace();
  auto& local_a=workspace.local_a;auto& local_b=workspace.local_b;
  CHECK(upload(*na,ar,ac,av,local_a,workspace.host));CHECK(upload(*nb,br,bc,bv,local_b,workspace.host));
  if(*nc)CHECK(upload(*nc,cr,cc,cv,workspace.local_c,workspace.host));
  CHECK(compute_product(rows,columns,n,row_offset,column_offset,na,local_a,nb,local_b,nc,workspace.local_c,
      alpha,beta,threshold,count,workspace,nullptr));
  Result* result=new Result;
  try {result->entries.resize(*count);}catch(...){delete result;throw;}
  auto error=cudaMemcpy(result->entries.data(),workspace.packed.p,size_t(*count)*sizeof(Coordinate),cudaMemcpyDeviceToHost);
  if(error!=cudaSuccess){delete result;return -1;}*output=result;
  return 0;
 }catch(const std::exception& e){std::fprintf(stderr,"Complex product: %s\n",e.what());return -1;}
}
extern "C" void ntpoly_complex_product_export_(void** handle,int* rows,int* columns,cuDoubleComplex* values) {
 auto* result=static_cast<Result*>(*handle);
 for(size_t i=0;i<result->entries.size();i++){rows[i]=result->entries[i].row;columns[i]=result->entries[i].column;values[i]=result->entries[i].value;}
 delete result;*handle=nullptr;
}

// Regression observability: scratch ownership and actual device growth only.
extern "C" int ntpoly_complex_workspace_count_(){return workspace_count;}
extern "C" uint64_t ntpoly_complex_allocation_count_(){return allocation_count;}

// Benchmark the actual shared-product implementation with a single input import
// and a single final export. Per-call times include panel communication, sort,
// accumulation, thresholding, sparse compaction, scalar nnz metadata and sync.
// No matrix payload crosses the CPU boundary inside the measured region.
extern "C" int ntpoly_complex_resident_benchmark_(int* rows,int* columns,int* n,
 int* row_offset,int* column_offset,int* nnz,const int* ir,const int* ic,const cuDoubleComplex* iv,
 double* threshold,int* repetitions,int* profile,double* seconds,double* phases,void** output,int* count) {
 *output=nullptr;*count=0;
 if(*repetitions<2)return -1;
 try {
  auto& workspace=get_workspace();
  CHECK(upload(*nnz,ir,ic,iv,workspace.local_a,workspace.host));
  CUDA(cudaDeviceSynchronize());
  int zero=0;double alpha=1,beta=0;
  for(int rep=0;rep<*repetitions;++rep) {
   if(*profile && rep==1)CUDA(cudaProfilerStart());
   auto start=std::chrono::steady_clock::now();
   CHECK(compute_product(rows,columns,n,row_offset,column_offset,nnz,workspace.local_a,
       nnz,workspace.local_a,&zero,workspace.local_c,&alpha,&beta,threshold,count,workspace,phases+4*rep));
   seconds[rep]=elapsed(start);
  }
  if(*profile)CUDA(cudaProfilerStop());
  auto result=std::make_unique<Result>();result->entries.resize(*count);
  CUDA(cudaMemcpy(result->entries.data(),workspace.packed.p,size_t(*count)*sizeof(Coordinate),cudaMemcpyDeviceToHost));
  *output=result.release();return 0;
 }catch(const std::exception& e){std::fprintf(stderr,"Complex resident benchmark: %s\n",e.what());return -1;}
}
