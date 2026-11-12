// Real FP64 Matrix_ps eigensolver boundary. No benchmark/file-format dependency.
#include <mpi.h>
#include <cuda_runtime.h>
#include <nccl.h>
#include <cusolverMp.h>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <string>
#include <vector>

namespace {
struct Coordinate { int row,column; double value; };
struct Entry { uint64_t offset; double value; };
struct Result { std::vector<Coordinate> entries; std::vector<double> eigenvalues; };
void check(int status,MPI_Comm comm,const char* operation) {
 if(status) { int rank;MPI_Comm_rank(comm,&rank);std::fprintf(stderr,"NTPoly cuSOLVERMp rank %d: %s failed (%d)\n",rank,operation,status);MPI_Abort(comm,1);std::abort(); }
}
#define CHECK(expr) check(int(expr),comm,#expr)
bool unanimous(bool value,MPI_Comm comm) { int a=value,b=0;MPI_Allreduce(&a,&b,1,MPI_INT,MPI_MIN,comm);return b!=0; }
struct Session {
 cudaStream_t stream=nullptr; ncclComm_t nccl=nullptr; cusolverMpHandle_t solver=nullptr;
 cusolverMpGrid_t grid=nullptr; cusolverMpMatrixDescriptor_t desc=nullptr;
 double *a=nullptr,*q=nullptr,*w=nullptr; void *work=nullptr; Entry* input=nullptr;
 int* info=nullptr;unsigned long long* count=nullptr;Coordinate* output=nullptr;
 ~Session() {
  if(stream)cudaStreamSynchronize(stream);
  for(void* p:{static_cast<void*>(a),static_cast<void*>(q),static_cast<void*>(w),work,
              static_cast<void*>(input),static_cast<void*>(info),static_cast<void*>(count),static_cast<void*>(output)})if(p)cudaFree(p);
  if(desc)cusolverMpDestroyMatrixDesc(desc);
  if(grid)cusolverMpDestroyGrid(grid);
  if(solver)cusolverMpDestroy(solver);
  if(nccl)ncclCommDestroy(nccl);
  if(stream)cudaStreamDestroy(stream);
 }
};
__global__ void fill(double* a,const Entry* entries,size_t count) {
 for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<count;i+=size_t(blockDim.x)*gridDim.x)
  a[entries[i].offset]=entries[i].value;
}
__global__ void compact(const double* q,int64_t rows,int64_t columns,int64_t ld,
 int block,int pr,int pc,int nr,int nc,int nvals,double threshold,unsigned long long* count,Coordinate* out) {
 unsigned long long total=0;
 for(int64_t x=int64_t(blockIdx.x)*blockDim.x+threadIdx.x;x<rows*columns;x+=int64_t(blockDim.x)*gridDim.x) {
  int64_t i=x%rows,j=x/rows;
  int r=int((i/block*nr+pr)*block+i%block),c=int((j/block*nc+pc)*block+j%block);
  double v=q[i+j*ld];
  if(c<nvals&&fabs(v)>threshold) {
   if(out) {auto index=atomicAdd(count,1ULL);out[index]={r+1,c+1,v};}else ++total;
  }
 }
 if(!out)atomicAdd(count,total);
}
// Capacity failures decline collectively. Runtime/numerical failures are errors,
// not a silent switch to a different answer after a failed collective.
int solve(MPI_Comm comm,int n,int nr,int nc,int block,int nvals,bool vectors,double threshold,
          const std::vector<Coordinate>& coordinates,Result& result) {
 int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
 MPI_Comm local;CHECK(MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local));
 int local_rank,local_size;MPI_Comm_rank(local,&local_rank);MPI_Comm_size(local,&local_size);
 const char* multi=std::getenv("YUKAWA_NCCL_ALLOW_MULTINODE");
 bool available=(local_size==size&&size<=4)||(multi&&std::string(multi)=="1");
 int devices=0;available=available&&cudaGetDeviceCount(&devices)==cudaSuccess&&devices>0;
 const char* env=std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");if(!env)env=std::getenv("SLURM_LOCALID");
 int device=devices?(env?std::atoi(env):local_rank)%devices:0;
 cudaDeviceProp properties{};
 available=available&&cudaGetDeviceProperties(&properties,device)==cudaSuccess;
 std::vector<unsigned char> identities(size_t(local_size)*16);
 CHECK(MPI_Allgather(properties.uuid.bytes,16,MPI_BYTE,identities.data(),16,MPI_BYTE,local));
 for(int i=0;i<local_size;i++)for(int j=i+1;j<local_size;j++)
  if(std::equal(identities.begin()+16*i,identities.begin()+16*(i+1),identities.begin()+16*j))available=false;
 MPI_Comm_free(&local);
 if(!unanimous(available,comm))return 1;
 CHECK(cudaSetDevice(device));
 int pr=rank/nc,pc=rank%nc;
 int64_t rows=cusolverMpNUMROC(n,block,pr,0,nr),cols=cusolverMpNUMROC(n,block,pc,0,nc),ld=std::max<int64_t>(1,rows);
 // Installed0.7.2 gives an incorrect spectrum for an empty-owner layout
 // (N129, block64,1x4), despite numerical info0. Decline collectively.
 if(!unanimous(rows>0&&cols>0&&(!vectors||uint64_t(rows)*uint64_t(cols)<=INT_MAX),comm))return 1;
 size_t cells=size_t(ld)*std::max<int64_t>(1,cols),bytes=cells*sizeof(double),free_mem,total_mem;
 CHECK(cudaMemGetInfo(&free_mem,&total_mem));
 size_t cap=std::numeric_limits<size_t>::max();const char* cap_env=std::getenv("NTPOLY_GPU_EIGEN_MEMORY_MB");
 if(cap_env&&*cap_env) {
  char* end=nullptr;unsigned long long mb=std::strtoull(cap_env,&end,10);
  check((*cap_env=='-'||*end||mb>SIZE_MAX/1048576)?1:0,comm,"invalid NTPOLY_GPU_EIGEN_MEMORY_MB");cap=size_t(mb)*1048576;
 }
 size_t fixed=2*bytes+size_t(n)*8+coordinates.size()*sizeof(Entry)+sizeof(int)+sizeof(unsigned long long);
 size_t export_bound=vectors?size_t(rows)*size_t(cols)*sizeof(Coordinate):0;
 if(!unanimous(fixed<=cap&&export_bound<=cap-fixed&&fixed+export_bound<=free_mem*0.8,comm))return 1;
 Session s;
 CHECK(cudaStreamCreate(&s.stream));ncclUniqueId id;if(rank==0)CHECK(ncclGetUniqueId(&id));
 CHECK(MPI_Bcast(&id,sizeof(id),MPI_BYTE,0,comm));CHECK(ncclCommInitRank(&s.nccl,size,id,rank));
 CHECK(cusolverMpCreate(&s.solver,device,s.stream));
 CHECK(cusolverMpCreateDeviceGrid(s.solver,&s.grid,s.nccl,nr,nc,CUSOLVERMP_GRID_MAPPING_ROW_MAJOR));
 CHECK(cusolverMpCreateMatrixDesc(&s.desc,s.grid,CUDA_R_64F,n,n,block,block,0,0,ld));
 bool allocated=cudaMalloc(&s.a,bytes)==cudaSuccess;
 allocated=(cudaMalloc(&s.q,bytes)==cudaSuccess)&&allocated;
 allocated=(cudaMalloc(&s.w,size_t(n)*8)==cudaSuccess)&&allocated;
 allocated=(cudaMalloc(&s.info,sizeof(int))==cudaSuccess)&&allocated;
 if(!unanimous(allocated,comm)) {cudaGetLastError();return 1;}
 char compz='V';size_t dw=0,hw=0;
 CHECK(cusolverMpSyevd_bufferSize(s.solver,&compz,CUBLAS_FILL_MODE_LOWER,n,s.a,1,1,s.desc,s.w,s.q,1,1,s.desc,CUDA_R_64F,&dw,&hw));
 CHECK(cudaMemGetInfo(&free_mem,&total_mem));
 if(!unanimous(dw<=cap-fixed-export_bound&&dw+export_bound<=free_mem*0.8,comm))return 1;
 std::vector<unsigned char> host;
 std::vector<Entry> packed;
 try {
  host.resize(std::max<size_t>(hw,1));packed.reserve(coordinates.size());
  for(auto e:coordinates) {
   int r=e.row-1,c=e.column-1;
   packed.push_back({uint64_t((r/block/nr)*block+r%block)+uint64_t((c/block/nc)*block+c%block)*ld,e.value});
  }
 } catch(const std::bad_alloc&) {allocated=false;}
 allocated=(cudaMalloc(&s.work,std::max<size_t>(dw,1))==cudaSuccess)&&allocated;
 allocated=(cudaMalloc(&s.input,std::max<size_t>(1,packed.size())*sizeof(Entry))==cudaSuccess)&&allocated;
 if(!unanimous(allocated,comm)){cudaGetLastError();return 1;}
 CHECK(cudaMemcpyAsync(s.input,packed.data(),packed.size()*sizeof(Entry),cudaMemcpyHostToDevice,s.stream));
 CHECK(cudaMemsetAsync(s.a,0,bytes,s.stream));fill<<<256,256,0,s.stream>>>(s.a,s.input,packed.size());CHECK(cudaGetLastError());
 CHECK(cusolverMpSyevd(s.solver,&compz,CUBLAS_FILL_MODE_LOWER,n,s.a,1,1,s.desc,s.w,s.q,1,1,s.desc,CUDA_R_64F,s.work,dw,host.data(),hw,s.info));
 CHECK(cudaStreamSynchronize(s.stream));int info;CHECK(cudaMemcpy(&info,s.info,sizeof(int),cudaMemcpyDeviceToHost));
 check(info,comm,"Syevd numerical info");
 if(rank==0){result.eigenvalues.resize(nvals);CHECK(cudaMemcpy(result.eigenvalues.data(),s.w,size_t(nvals)*8,cudaMemcpyDeviceToHost));}
 if(vectors) {
  CHECK(cudaMalloc(&s.count,sizeof(unsigned long long)));CHECK(cudaMemsetAsync(s.count,0,sizeof(unsigned long long),s.stream));
  compact<<<256,256,0,s.stream>>>(s.q,rows,cols,ld,block,pr,pc,nr,nc,nvals,threshold,s.count,nullptr);CHECK(cudaGetLastError());
  CHECK(cudaStreamSynchronize(s.stream));unsigned long long count;CHECK(cudaMemcpy(&count,s.count,sizeof(count),cudaMemcpyDeviceToHost));
  check(count>INT_MAX?1:0,comm,"eigenvector output exceeds Matrix_ps local index capacity");
  result.entries.resize(size_t(count));CHECK(cudaMalloc(&s.output,std::max<size_t>(count,1)*sizeof(Coordinate)));
  CHECK(cudaMemsetAsync(s.count,0,sizeof(unsigned long long),s.stream));
  compact<<<256,256,0,s.stream>>>(s.q,rows,cols,ld,block,pr,pc,nr,nc,nvals,threshold,s.count,s.output);CHECK(cudaGetLastError());
  CHECK(cudaStreamSynchronize(s.stream));CHECK(cudaMemcpy(result.entries.data(),s.output,size_t(count)*sizeof(Coordinate),cudaMemcpyDeviceToHost));
 }
 if(rank==0)std::printf("NTPoly cuSOLVERMp eigen: n=%d active=%d grid=%dx%d block=%d nvals=%d vectors=%d\n",n,size,nr,nc,block,nvals,int(vectors));
 return 0;
}
}

extern "C" int ntpoly_cusolvermp_eigh_(int* fcomm,int* n,int* nr,int* nc,int* count,
 const int* rows,const int* columns,const double* values,int* nvals,int* vectors,double* threshold,int* root,void** out) {
 MPI_Comm comm=MPI_Comm_f2c(MPI_Fint(*fcomm));int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
 *out=nullptr;
 try {
  if(uint64_t(*n)*uint64_t(*n)>SIZE_MAX/64)return 1;
  int block=*n<1024?64:512,pr=1,pc=1;
  if(!*root) {
   for(int candidate=1;candidate<=size/candidate;candidate++)if(size%candidate==0)pr=candidate;
   pc=size/pr;
  }
  MPI_Datatype coordinate_type;CHECK(MPI_Type_contiguous(sizeof(Coordinate),MPI_BYTE,&coordinate_type));CHECK(MPI_Type_commit(&coordinate_type));
  std::vector<int> counts(size,0),recv_counts(size),offsets(size),recv_offsets(size),cursor;
  std::vector<int64_t> planned(size,0);
  std::vector<Coordinate> input;
  // Match NTPoly's LAPACK upper-triangle convention. Materialize its symmetric
  // counterpart as well for the distributed implementation.
  auto visit=[&](auto append) {
   for(int i=0;i<*count;i++) {
    check(rows[i]<1||columns[i]<1||!std::isfinite(values[i]),comm,"invalid eigensolver input");
    if(rows[i]>*n||columns[i]>*n||rows[i]>columns[i])continue;
    append(Coordinate{rows[i],columns[i],values[i]});
    if(rows[i]!=columns[i])append(Coordinate{columns[i],rows[i],values[i]});
   }
  };
  auto owner=[&](Coordinate e){return *root?0:((e.row-1)/block%pr)*pc+(e.column-1)/block%pc;};
  visit([&](Coordinate e){planned[owner(e)]++;});
  bool counts_safe=true;
  for(int i=0;i<size;i++){if(planned[i]>INT_MAX)counts_safe=false;else counts[i]=int(planned[i]);}
  if(!unanimous(counts_safe,comm)){MPI_Type_free(&coordinate_type);return 1;}
  CHECK(MPI_Alltoall(counts.data(),1,MPI_INT,recv_counts.data(),1,MPI_INT,comm));
  int64_t total=0,in_total=0;bool safe=true;
  for(int i=0;i<size;i++){offsets[i]=int(total);recv_offsets[i]=int(in_total);total+=counts[i];in_total+=recv_counts[i];if(total>INT_MAX||in_total>INT_MAX)safe=false;}
  if(!unanimous(safe,comm)){MPI_Type_free(&coordinate_type);return 1;}
  input.resize(size_t(total));cursor=offsets;
  visit([&](Coordinate e){input[cursor[owner(e)]++]=e;});
  std::vector<Coordinate> received(size_t(in_total),Coordinate{});
  CHECK(MPI_Alltoallv(input.data(),counts.data(),offsets.data(),coordinate_type,received.data(),recv_counts.data(),recv_offsets.data(),coordinate_type,comm));
  CHECK(MPI_Type_free(&coordinate_type));
  Result* result=new Result;int status=0;
  if(!*root||rank==0)status=solve(*root?MPI_COMM_SELF:comm,*n,pr,pc,block,*nvals,*vectors!=0,*threshold,received,*result);
  if(*root)CHECK(MPI_Bcast(&status,1,MPI_INT,0,comm));
  if(status){delete result;return status;}
  *out=result;return 0;
 } catch(const std::bad_alloc&) {check(1,comm,"host eigensolver allocation");return -1;}
}
extern "C" void ntpoly_cusolvermp_result_(void** handle,int* count,const double** eigenvalues,int* eigen_count) {
 auto* result=static_cast<Result*>(*handle);*count=int(result->entries.size());*eigen_count=int(result->eigenvalues.size());
 *eigenvalues=result->eigenvalues.data();
}
extern "C" void ntpoly_cusolvermp_coordinates_(void** handle,int* rows,int* columns,double* values) {
 const auto& entries=static_cast<Result*>(*handle)->entries;
 for(size_t i=0;i<entries.size();i++){rows[i]=entries[i].row;columns[i]=entries[i].column;values[i]=entries[i].value;}
}
extern "C" void ntpoly_cusolvermp_free_(void** handle){delete static_cast<Result*>(*handle);*handle=nullptr;}
