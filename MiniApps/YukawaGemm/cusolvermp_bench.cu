// Distributed FP64 diagonalization and occupied density benchmark.
// Standalone build/run instructions and timing scopes: CUSOLVERMP.md.
#include <mpi.h>
#include <cuda_runtime.h>
#include <nccl.h>
#include <cusolverMp.h>
#include <cublasmp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <vector>

static int rank_id;
#define CHECK(call) do { int status = int(call); if(status != 0) { \
 fprintf(stderr,"rank %d: %s returned %d at %d\n",rank_id,#call,status,__LINE__); \
 MPI_Abort(MPI_COMM_WORLD,1); } } while(0)
static void require(bool ok, const char* message) {
 if(!ok) { fprintf(stderr,"rank %d: %s\n",rank_id,message); MPI_Abort(MPI_COMM_WORLD,2); }
}
static double elapsed(double start) {
 CHECK(cudaDeviceSynchronize()); double t=MPI_Wtime()-start, maximum;
 CHECK(MPI_Allreduce(&t,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD)); return maximum;
}
struct Entry { uint64_t offset; double value; };
struct Coordinate { int32_t col,row; double value; };
__global__ void compact(const double* a,int64_t rows,int64_t cols,int64_t ld,
 int block,int pr,int pc,int nr,int nc,unsigned long long* count,Coordinate* output) {
 unsigned long long total=0;
 for(int64_t x=int64_t(blockIdx.x)*blockDim.x+threadIdx.x;x<rows*cols;x+=int64_t(blockDim.x)*gridDim.x) {
  int64_t i=x%rows,j=x/rows;double v=a[i+j*ld];
  if(fabs(v)>1e-4) {
   if(output) {auto index=atomicAdd(count,1ULL);output[index]={int32_t((j/block*nc+pc)*block+j%block+1),int32_t((i/block*nr+pr)*block+i%block+1),v};}
   else ++total;
  }
 }
 if(!output)atomicAdd(count,total);
}
__global__ void fill(double* a,const Entry* entries,size_t count) {
 for(size_t i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=size_t(blockDim.x)*gridDim.x)
  a[entries[i].offset]=entries[i].value;
}
// Atomics operate on block sums to avoid one atomic per matrix element.
__global__ void metrics(const double* a,int64_t rows,int64_t cols,int64_t ld,
 int block,int pr,int pc,int nr,int nc,double* result) {
 __shared__ double sums[3][256]; double tr=0,all=0,sq=0;
 for(int64_t x=int64_t(blockIdx.x)*blockDim.x+threadIdx.x;x<rows*cols;x+=int64_t(blockDim.x)*gridDim.x) {
  int64_t i=x%rows,j=x/rows;
  int64_t gi=(i/block*nr+pr)*block+i%block;
  int64_t gj=(j/block*nc+pc)*block+j%block;
  double v=a[i+j*ld]; all+=v; sq+=v*v; if(gi==gj)tr+=v;
 }
 sums[0][threadIdx.x]=tr;sums[1][threadIdx.x]=all;sums[2][threadIdx.x]=sq;
 __syncthreads();
 for(int stride=128;stride;stride/=2) { if(threadIdx.x<stride)
  for(int k=0;k<3;k++)sums[k][threadIdx.x]+=sums[k][threadIdx.x+stride]; __syncthreads(); }
 if(threadIdx.x==0)for(int k=0;k<3;k++)atomicAdd(result+k,sums[k][0]);
}
__global__ void energy(const double* a,const Entry* e,size_t count,double* result) {
 __shared__ double sums[256]; double s=0;
 for(size_t i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=size_t(blockDim.x)*gridDim.x)
  s+=a[e[i].offset]*e[i].value;
 sums[threadIdx.x]=s;__syncthreads();
 for(int stride=128;stride;stride/=2) {if(threadIdx.x<stride)sums[threadIdx.x]+=sums[threadIdx.x+stride];__syncthreads();}
 if(threadIdx.x==0)atomicAdd(result+3,sums[0]);
}
int main(int argc,char**argv) {
 CHECK(MPI_Init(&argc,&argv)); int size;CHECK(MPI_Comm_rank(MPI_COMM_WORLD,&rank_id));CHECK(MPI_Comm_size(MPI_COMM_WORLD,&size));
 require(argc==6||argc==7,"usage: cusolvermp_bench input.bin prows pcols nelec repetitions [block=256]");
 int nr=atoi(argv[2]),nc=atoi(argv[3]),occ=atoi(argv[4]),reps=atoi(argv[5]),block=argc==7?atoi(argv[6]):256;
 require(nr>0&&nc>0&&nr*nc==size&&reps>0&&block>0,"invalid grid, block or repetitions");
 double start=MPI_Wtime(); MPI_Comm local;CHECK(MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank_id,MPI_INFO_NULL,&local));
 int lr,devices;CHECK(MPI_Comm_rank(local,&lr));CHECK(cudaGetDeviceCount(&devices));require(lr<devices,"one MPI rank per GPU required");CHECK(cudaSetDevice(lr));
 cudaStream_t stream;CHECK(cudaStreamCreate(&stream));ncclUniqueId id;
 if(!rank_id)CHECK(ncclGetUniqueId(&id));CHECK(MPI_Bcast(&id,sizeof(id),MPI_BYTE,0,MPI_COMM_WORLD));
 ncclComm_t comm;CHECK(ncclCommInitRank(&comm,size,id,rank_id));
 cusolverMpHandle_t solver;CHECK(cusolverMpCreate(&solver,lr,stream));
 cusolverMpGrid_t grid;CHECK(cusolverMpCreateDeviceGrid(solver,&grid,comm,nr,nc,CUSOLVERMP_GRID_MAPPING_ROW_MAJOR));
 cublasMpHandle_t blas;CHECK(cublasMpCreate(&blas,stream));cublasMpGrid_t bg;
 CHECK(cublasMpGridCreate(nr,nc,CUBLASMP_GRID_LAYOUT_ROW_MAJOR,comm,&bg));
 int version;CHECK(cusolverMpGetVersion(solver,&version));
 std::ifstream f(argv[1],std::ios::binary);int32_t n,m,complex;int64_t nnz;
 f.read((char*)&n,4);f.read((char*)&m,4);f.read((char*)&complex,4);f.read((char*)&nnz,8);
 require(bool(f)&&n==m&&n>0&&complex==0&&nnz>=0&&occ>0&&occ<n,"invalid real symmetric input header/electron count");
 int pr=rank_id/nc,pc=rank_id%nc;
 int64_t rows=cusolverMpNUMROC(n,block,pr,0,nr),cols=cusolverMpNUMROC(n,block,pc,0,nc),ld=std::max<int64_t>(1,rows);
 size_t bytes=ld*std::max<int64_t>(1,cols)*sizeof(double); std::vector<Entry> entries;
 // Input is a full symmetric coordinate matrix, with one-based column,row,value records.
 for(int64_t k=0;k<nnz;k++) { int32_t c,r;double v;f.read((char*)&c,4);f.read((char*)&r,4);f.read((char*)&v,8);
  require(bool(f)&&r>=1&&r<=n&&c>=1&&c<=n&&std::isfinite(v),"invalid coordinate record");--r;--c;
  if((r/block)%nr==pr&&(c/block)%nc==pc)entries.push_back({uint64_t((r/block/nr)*block+r%block)+uint64_t((c/block/nc)*block+c%block)*ld,v});
 }
 double *a,*q,*eval,*stats;Entry* de; int* info;
 size_t free_mem,total_mem;CHECK(cudaMemGetInfo(&free_mem,&total_mem));
 require(2*bytes+size_t(n)*8<free_mem*0.8,"dense A and Q exceed 80% free GPU memory");
 CHECK(cudaMalloc(&a,bytes));CHECK(cudaMalloc(&q,bytes));CHECK(cudaMalloc(&eval,size_t(n)*8));
 CHECK(cudaMalloc(&de,std::max<size_t>(1,entries.size())*sizeof(Entry)));CHECK(cudaMemcpy(de,entries.data(),entries.size()*sizeof(Entry),cudaMemcpyHostToDevice));
 CHECK(cudaMalloc(&stats,4*sizeof(double)));CHECK(cudaMalloc(&info,sizeof(int)));
 cusolverMpMatrixDescriptor_t desc;CHECK(cusolverMpCreateMatrixDesc(&desc,grid,CUDA_R_64F,n,n,block,block,0,0,ld));
 cublasMpMatrixDescriptor_t bd;CHECK(cublasMpMatrixDescriptorCreate(n,n,block,block,0,0,ld,CUDA_R_64F,bg,&bd));
 char compz='V';size_t sd,sh,gd,gh;double one=1,zero=0;
 CHECK(cusolverMpSyevd_bufferSize(solver,&compz,CUBLAS_FILL_MODE_LOWER,n,a,1,1,desc,eval,q,1,1,desc,CUDA_R_64F,&sd,&sh));
 CHECK(cublasMpGemm_bufferSize(blas,CUBLAS_OP_N,CUBLAS_OP_T,n,n,occ,&one,q,1,1,bd,q,1,1,bd,&zero,a,1,1,bd,CUBLAS_COMPUTE_64F,&gd,&gh));
 size_t dw=std::max(sd,gd),hw=std::max(sh,gh);CHECK(cudaMemGetInfo(&free_mem,&total_mem));
 require(dw<free_mem*0.8,"workspace exceeds 80% remaining GPU memory");
 void* work;CHECK(cudaMalloc(&work,std::max<size_t>(dw,1))); std::vector<unsigned char> host(std::max<size_t>(hw,1));
 double setup=elapsed(start);size_t peak=2*bytes+dw+size_t(n)*8+entries.size()*sizeof(Entry);
 size_t peakmax;CHECK(MPI_Allreduce(&peak,&peakmax,1,MPI_UINT64_T,MPI_MAX,MPI_COMM_WORLD));
 if(!rank_id)printf("config n=%d ranks=%d grid=%dx%d block=%d solver_version=%d setup=%.9g max_rank_explicit_device_bytes=%zu solver_work=%zu gemm_work=%zu host_work=%zu\n",n,size,nr,nc,block,version,setup,peakmax,sd,gd,hw);
 for(int rep=0;rep<reps;rep++) {
  CHECK(MPI_Barrier(MPI_COMM_WORLD));start=MPI_Wtime();
  CHECK(cudaMemsetAsync(a,0,bytes,stream));fill<<<256,256,0,stream>>>(a,de,entries.size());CHECK(cudaGetLastError());double import=elapsed(start);
  CHECK(MPI_Barrier(MPI_COMM_WORLD));start=MPI_Wtime();
  CHECK(cusolverMpSyevd(solver,&compz,CUBLAS_FILL_MODE_LOWER,n,a,1,1,desc,eval,q,1,1,desc,CUDA_R_64F,work,sd,host.data(),sh,info));
  double evd=elapsed(start);int result;CHECK(cudaMemcpy(&result,info,sizeof(int),cudaMemcpyDeviceToHost));require(result==0,"Syevd numerical failure");
  CHECK(MPI_Barrier(MPI_COMM_WORLD));start=MPI_Wtime();
  CHECK(cublasMpGemm(blas,CUBLAS_OP_N,CUBLAS_OP_T,n,n,occ,&one,q,1,1,bd,q,1,1,bd,&zero,a,1,1,bd,CUBLAS_COMPUTE_64F,work,gd,host.data(),gh));
  double density=elapsed(start);
  start=MPI_Wtime();CHECK(cudaMemsetAsync(stats,0,4*sizeof(double),stream));
  metrics<<<256,256,0,stream>>>(a,rows,cols,ld,block,pr,pc,nr,nc,stats);energy<<<256,256,0,stream>>>(a,de,entries.size(),stats);CHECK(cudaGetLastError());CHECK(cudaStreamSynchronize(stream));
  double s[4],global[4];CHECK(cudaMemcpy(s,stats,sizeof(s),cudaMemcpyDeviceToHost));CHECK(MPI_Allreduce(s,global,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD));
  std::vector<double> e(n);CHECK(cudaMemcpy(e.data(),eval,size_t(n)*8,cudaMemcpyDeviceToHost));double esum=0;for(int i=0;i<occ;i++)esum+=e[i];
  double validation=elapsed(start);
  require(std::abs(global[0]-occ)<1e-7*occ&&std::abs(global[2]-occ)<1e-7*occ&&std::abs(global[3]-esum)<1e-7*std::max(1.,std::abs(esum)),"density trace/norm/energy validation failed");
  // Exhaustive CPU checks for small validation inputs, outside solve timing.
  if(n<=1024) {
   std::vector<double> local_q(ld*std::max<int64_t>(1,cols)),full_q(size_t(n)*n,0),full_h(size_t(n)*n,0),full_d(size_t(n)*n,0),local_d(local_q.size());
   CHECK(cudaMemcpy(local_q.data(),q,bytes,cudaMemcpyDeviceToHost));CHECK(cudaMemcpy(local_d.data(),a,bytes,cudaMemcpyDeviceToHost));
   for(int64_t j=0;j<cols;j++)for(int64_t i=0;i<rows;i++) {
    int64_t gi=(i/block*nr+pr)*block+i%block,gj=(j/block*nc+pc)*block+j%block;
    full_q[gi+gj*n]=local_q[i+j*ld];full_d[gi+gj*n]=local_d[i+j*ld];
   }
   for(auto entry:entries) {
    int64_t i=entry.offset%ld,j=entry.offset/ld;
    int64_t gi=(i/block*nr+pr)*block+i%block,gj=(j/block*nc+pc)*block+j%block;
    full_h[gi+gj*n]=entry.value;
   }
   for(auto* matrix:{&full_q,&full_h,&full_d})CHECK(MPI_Allreduce(MPI_IN_PLACE,matrix->data(),n*n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD));
   if(!rank_id) {
    double residual=0,orthogonal=0,projector=0,symmetry=0,hnorm=0;
    for(double v:full_h)hnorm+=v*v;
    for(int j=0;j<n;j++)for(int i=0;i<n;i++) {
     double hq=0,qq=0,p=0;
     for(int k=0;k<n;k++) {hq+=full_h[i+k*n]*full_q[k+j*n];qq+=full_q[k+i*n]*full_q[k+j*n];if(k<occ)p+=full_q[i+k*n]*full_q[j+k*n];}
     double r=hq-full_q[i+j*n]*e[j];residual+=r*r;double o=qq-(i==j);orthogonal+=o*o;
     projector=std::max(projector,std::abs(p-full_d[i+j*n]));symmetry=std::max(symmetry,std::abs(full_h[i+j*n]-full_h[j+i*n]));
    }
    residual=std::sqrt(residual/std::max(1.,hnorm));orthogonal=std::sqrt(orthogonal);
    printf("small_check residual=%.9g orthogonality=%.9g density_max_error=%.9g input_asymmetry=%.9g\n",residual,orthogonal,projector,symmetry);
    require(residual<1e-10&&orthogonal<1e-10&&projector<1e-10&&symmetry<1e-14,"exhaustive small matrix check failed");
   }
  }
  CHECK(MPI_Barrier(MPI_COMM_WORLD));start=MPI_Wtime();
  unsigned long long* count;CHECK(cudaMalloc(&count,sizeof(*count)));CHECK(cudaMemsetAsync(count,0,sizeof(*count),stream));
  compact<<<256,256,0,stream>>>(a,rows,cols,ld,block,pr,pc,nr,nc,count,nullptr);CHECK(cudaStreamSynchronize(stream));
  unsigned long long kept,global_kept;CHECK(cudaMemcpy(&kept,count,sizeof(kept),cudaMemcpyDeviceToHost));
  CHECK(cudaMemGetInfo(&free_mem,&total_mem));require(kept*sizeof(Coordinate)<free_mem*0.8,"sparse density output exceeds memory budget");
  Coordinate* output;CHECK(cudaMalloc(&output,std::max<size_t>(1,kept)*sizeof(Coordinate)));CHECK(cudaMemsetAsync(count,0,sizeof(*count),stream));
  compact<<<256,256,0,stream>>>(a,rows,cols,ld,block,pr,pc,nr,nc,count,output);CHECK(cudaGetLastError());CHECK(cudaStreamSynchronize(stream));
  std::vector<Coordinate> host_density(kept);CHECK(cudaMemcpy(host_density.data(),output,kept*sizeof(Coordinate),cudaMemcpyDeviceToHost));
  CHECK(cudaFree(output));CHECK(cudaFree(count));double export_seconds=elapsed(start);
  CHECK(MPI_Allreduce(&kept,&global_kept,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD));
  if(!rank_id) {printf("result rep=%d import=%.9g evd=%.9g density=%.9g solve=%.9g export=%.9g complete=%.9g validation=%.9g output_nnz=%llu energy=%.17g trace=%.17g sum=%.17g norm2=%.17g eig_energy=%.17g gap=%.17g info=%d\n",rep,import,evd,density,evd+density,export_seconds,import+evd+density+export_seconds,validation,global_kept,global[3],global[0],global[1],global[2],esum,e[occ]-e[occ-1],result);fflush(stdout);}
 }
 CHECK(cudaFree(work));CHECK(cudaFree(a));CHECK(cudaFree(q));CHECK(cudaFree(eval));CHECK(cudaFree(de));CHECK(cudaFree(stats));CHECK(cudaFree(info));
 CHECK(cublasMpMatrixDescriptorDestroy(bd));CHECK(cusolverMpDestroyMatrixDesc(desc));CHECK(cublasMpGridDestroy(bg));CHECK(cusolverMpDestroyGrid(grid));CHECK(cublasMpDestroy(blas));CHECK(cusolverMpDestroy(solver));CHECK(ncclCommDestroy(comm));CHECK(cudaStreamDestroy(stream));CHECK(MPI_Comm_free(&local));CHECK(MPI_Finalize());
}
