! Experimental owned GPU matrices for the library TRS2 iteration.
! Host Matrix_ps remains authoritative outside the explicitly bounded GPU region.
MODULE PSMatrixGPUBackendModule
  USE ISO_C_BINDING
  USE DataTypesModule, ONLY: NTREAL, MPINTLONG, MPINTREAL
  USE ISO_FORTRAN_ENV, ONLY: INT64, REAL64
  USE GPUPlanningModule, ONLY: EstimateGPUPlan, GPUMemoryBudget
  USE NTMPIModule
  USE PSMatrixModule, ONLY: Matrix_ps, ConstructEmptyMatrix, MergeMatrixLocalBlocks, SplitMatrixToLocalBlocks, &
    GatherMatrixToProcess, FillMatrixFromTripletList
  USE SMatrixModule, ONLY: Matrix_lsr, DestructMatrix, MatrixToTripletList
  USE TripletListModule, ONLY: TripletList_r, ConstructTripletList, DestructTripletList
  USE SolverCSR, ONLY: TransposeCSR
  USE SolverNCCL, ONLY: InitSolverNCCL, WarmSolverNCCL, FinishSolverNCCL, SolverScalarSum, nccl_active
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GPUMatrix_t, GPURequested, SelectGPUBackend, CheckGPUMemory, BeginGPUBackend, EndGPUBackend
  PUBLIC :: UploadGPURegion, DownloadGPURegion
  PUBLIC :: UploadGPU, DownloadGPU, CopyGPU, ScaleGPU, IncrementGPU, TraceGPU, DotGPU, SquareGPU, MultiplyGPU
  PUBLIC :: TryGPUSimilarityTransform, ReleaseGPU, BoundsGPURegion
  TYPE :: GPUMatrix_t
    PRIVATE
    INTEGER(C_INT) :: handle=0, rows=0, columns=0, dimension=0, diagonal=0
    INTEGER(C_INT) :: process_rows=0, process_columns=0
  CONTAINS
    FINAL :: DestroyGPU
    PROCEDURE, PRIVATE :: AssignGPU
    GENERIC :: ASSIGNMENT(=) => AssignGPU
  END TYPE
  INTERFACE
    INTEGER(C_INT) FUNCTION import_csc(slot,rows,columns,nnz,cp,ri,values) BIND(C,NAME='gpu_ops_import_csc_')
      IMPORT C_INT,C_PTR
      INTEGER(C_INT) :: slot,rows,columns,nnz
      TYPE(C_PTR), VALUE :: cp,ri,values
    END FUNCTION
    INTEGER(C_INT) FUNCTION bounds_slot(slot,columns,diagonal,lower,upper) BIND(C,NAME='gpu_ops_bounds_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT) :: slot,columns,diagonal
      REAL(C_DOUBLE) :: lower(*),upper(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION create_context(context) BIND(C,NAME='solver_nccl_context_create_')
      IMPORT C_PTR,C_INT
      TYPE(C_PTR) :: context
    END FUNCTION
    INTEGER(C_INT) FUNCTION select_context(context) BIND(C,NAME='solver_nccl_context_select_')
      IMPORT C_PTR,C_INT
      TYPE(C_PTR), VALUE :: context
    END FUNCTION
    INTEGER(C_INT) FUNCTION suspend_context() BIND(C,NAME='solver_nccl_context_suspend_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT) FUNCTION device_info(device,count,identity) BIND(C,NAME='gpu_device_info_')
      IMPORT C_INT,C_INT8_T
      INTEGER(C_INT) :: device,count
      INTEGER(C_INT8_T) :: identity(16)
    END FUNCTION
    INTEGER(C_INT) FUNCTION memory_info(free_bytes,total_bytes) BIND(C,NAME='gpu_memory_info_')
      IMPORT C_INT,C_INT64_T
      INTEGER(C_INT64_T) :: free_bytes,total_bytes
    END FUNCTION
    SUBROUTINE init_gpu(status,device) BIND(C,NAME='gpu_ops_init_')
      IMPORT C_INT
      INTEGER(C_INT) :: status,device
    END SUBROUTINE
    SUBROUTINE create_slot(slot) BIND(C,NAME='gpu_ops_create_')
      IMPORT C_INT
      INTEGER(C_INT) :: slot
    END SUBROUTINE
    SUBROUTINE destroy_slot(slot) BIND(C,NAME='gpu_ops_destroy_')
      IMPORT C_INT
      INTEGER(C_INT) :: slot
    END SUBROUTINE
    SUBROUTINE allocate_slot(slot,rows,nnz,rp,ci,values,status) BIND(C,NAME='gpu_ops_alloc_')
      IMPORT C_INT,C_PTR
      INTEGER(C_INT) :: slot,rows,nnz,status
      TYPE(C_PTR) :: rp,ci,values
    END SUBROUTINE
    SUBROUTINE get_slot(slot,rp,ci,values,rows,nnz) BIND(C,NAME='gpu_ops_get_')
      IMPORT C_INT,C_PTR
      INTEGER(C_INT) :: slot,rows,nnz
      TYPE(C_PTR) :: rp,ci,values
    END SUBROUTINE
    SUBROUTINE copy_slot(dst,src) BIND(C,NAME='gpu_ops_copy_')
      IMPORT C_INT
      INTEGER(C_INT) :: dst,src
    END SUBROUTINE
    SUBROUTINE scale_slot(slot,alpha) BIND(C,NAME='gpu_ops_scale_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT) :: slot
      REAL(C_DOUBLE) :: alpha
    END SUBROUTINE
    SUBROUTINE increment_slot(dst,src,alpha,threshold,nnz) BIND(C,NAME='gpu_ops_increment_threshold_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT) :: dst,src,nnz
      REAL(C_DOUBLE) :: alpha,threshold
    END SUBROUTINE
    SUBROUTINE trace_slot(slot,diagonal,value) BIND(C,NAME='gpu_ops_trace_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT) :: slot,diagonal
      REAL(C_DOUBLE) :: value
    END SUBROUTINE
    SUBROUTINE dot_slot(a,b,value) BIND(C,NAME='gpu_ops_dot_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT) :: a,b
      REAL(C_DOUBLE) :: value
    END SUBROUTINE
    INTEGER(C_INT) FUNCTION resident_square(input,output,m,n,k,pr,pc,threshold,phases) &
        BIND(C,NAME='gpu_resident_square_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT), VALUE :: input,output,m,n,k,pr,pc
      REAL(C_DOUBLE), VALUE :: threshold
      REAL(C_DOUBLE) :: phases(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION resident_multiply(left,right,output,m,n,k,pr,pc,threshold,phases) &
        BIND(C,NAME='gpu_resident_multiply_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT), VALUE :: left,right,output,m,n,k,pr,pc
      REAL(C_DOUBLE), VALUE :: threshold
      REAL(C_DOUBLE) :: phases(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION free_panels() BIND(C,NAME='gpu_resident_free_')
      IMPORT C_INT
    END FUNCTION
  END INTERFACE
CONTAINS
  SUBROUTINE FailGPU(message)
    CHARACTER(*), INTENT(IN) :: message
    INTEGER :: ierr
    PRINT *, 'NTPoly GPU backend: ',message
    CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    ERROR STOP
  END SUBROUTINE

  LOGICAL FUNCTION GPURequested(matrix) RESULT(requested)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    CHARACTER(16) :: option
    INTEGER :: flags(2),global_flags(2),ierr
    CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_TRS2',option)
    flags(1)=0
    IF(TRIM(option)=='1') flags(1)=1
    flags(2)=-flags(1)
    CALL MPI_Allreduce(flags,global_flags,2,MPI_INTEGER,MPI_MAX,matrix%process_grid%global_comm,ierr)
    IF(global_flags(1)/=-global_flags(2)) CALL FailGPU('NTPOLY_GPU_TRS2 differs between ranks')
    requested=global_flags(1)==1
    IF(.NOT.requested) RETURN
    IF(matrix%is_complex .OR. matrix%process_grid%num_process_slices/=1) &
      CALL FailGPU('TRS2 GPU backend currently requires real matrices and one slice')
  END FUNCTION

  ! Selection happens after CPU WH/X0 setup, using the actual working pattern.
  SUBROUTINE SelectGPUBackend(wh,x,enabled,root_only,expected_iterations,extra_fixed_nnz,initial_shift)
    TYPE(Matrix_ps), INTENT(IN) :: wh,x
    LOGICAL, INTENT(OUT) :: enabled,root_only
    INTEGER, INTENT(IN) :: expected_iterations
    INTEGER(INT64), OPTIONAL, INTENT(IN) :: extra_fixed_nnz
    LOGICAL, OPTIONAL, INTENT(IN) :: initial_shift
    CHARACTER(64) :: policy,cap_text,multinode
    INTEGER :: code,flags(2),global_flags(2),ierr,pr,pc,comm,rank,iterations,io,j,k
    INTEGER :: local_comm,local_size,valid,all_valid
    INTEGER(C_INT8_T), ALLOCATABLE :: devices(:,:)
    INTEGER(C_INT8_T) :: identity(16)
    INTEGER(C_INT) :: device,count,status
    INTEGER(INT64) :: nnz_wh,nnz_x,root_peak,all_peak,root_budget,all_budget,cap,local_cap
    INTEGER(INT64) :: cap_mb,local_wh,local_x,free_bytes,total_bytes
    REAL(REAL64) :: root_seconds,all_seconds
    LOGICAL :: root_index,all_index,root_fit,all_fit
    comm=wh%process_grid%global_comm; rank=wh%process_grid%global_rank
    enabled=.FALSE.; root_only=.FALSE.
    CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_RANK_POLICY',policy)
    SELECT CASE(TRIM(policy))
    CASE('','auto'); code=0
    CASE('all'); code=1
    CASE('root'); code=2
    CASE DEFAULT; CALL FailGPU('rank policy must be auto, all or root')
    END SELECT
    flags=[code,-code]
    CALL MPI_Allreduce(flags,global_flags,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    IF(global_flags(1)/=-global_flags(2)) CALL FailGPU('GPU rank policies differ between ranks')
    ! One application rank needs no gather/redistribution in automatic mode.
    IF(code==0 .AND. wh%process_grid%total_processors==1) code=1
    CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_MEMORY_MB',cap_text)
    local_cap=HUGE(local_cap)
    IF(LEN_TRIM(cap_text)>0) THEN
      READ(cap_text,*,IOSTAT=io) cap_mb
      IF(io/=0) CALL FailGPU('NTPOLY_GPU_MEMORY_MB must be a nonnegative integer MiB count')
      IF(cap_mb<0 .OR. cap_mb>HUGE(cap_mb)/1048576_INT64) CALL FailGPU('GPU memory cap is out of range')
      local_cap=cap_mb*1048576_INT64
    END IF
    CALL MPI_Allreduce(local_cap,cap,1,MPINTLONG,MPI_MIN,comm,ierr)
    IF(cap==HUGE(cap)) cap=-1_INT64
    ! Count local sub-blocks without merging/copying their payloads.
    local_wh=0; local_x=0
    DO j=1,SIZE(wh%local_data_r,2)
      DO k=1,SIZE(wh%local_data_r,1)
        local_wh=local_wh+SIZE(wh%local_data_r(k,j)%values,KIND=INT64)
        local_x=local_x+SIZE(x%local_data_r(k,j)%values,KIND=INT64)
      END DO
    END DO
    CALL MPI_Allreduce(local_wh,nnz_wh,1,MPINTLONG,MPI_SUM,comm,ierr)
    CALL MPI_Allreduce(local_x,nnz_x,1,MPINTLONG,MPI_SUM,comm,ierr)
    IF(PRESENT(initial_shift)) THEN
      IF(initial_shift) THEN
        ! X0 adds at most one diagonal per actual row. Use that upper bound
        ! for ranking without materializing X0 on the CPU first.
        nnz_x=MIN(nnz_x+INT(wh%actual_matrix_dimension,INT64), &
          INT(wh%logical_matrix_dimension,INT64)*INT(wh%logical_matrix_dimension,INT64))
      END IF
    END IF
    IF(PRESENT(extra_fixed_nnz)) THEN
      IF(extra_fixed_nnz<0 .OR. nnz_wh>HUGE(nnz_wh)-extra_fixed_nnz) CALL FailGPU('fixed input count overflow')
      nnz_wh=nnz_wh+extra_fixed_nnz
    END IF
    pr=wh%process_grid%num_process_rows; pc=wh%process_grid%num_process_columns
    iterations=MAX(1,MIN(24,expected_iterations))
    CALL EstimateGPUPlan(wh%actual_matrix_dimension,1,1,nnz_wh,nnz_x,iterations, &
      root_peak,root_seconds,root_index)
    CALL EstimateGPUPlan(wh%logical_matrix_dimension,pr,pc,nnz_wh,nnz_x,iterations, &
      all_peak,all_seconds,all_index)
    root_budget=0; all_budget=0; root_fit=.FALSE.; all_fit=.FALSE.
    IF(code/=1 .AND. root_index .AND. (cap<0 .OR. root_peak<=cap)) THEN
      IF(rank==0) THEN
        status=memory_info(free_bytes,total_bytes)
        IF(status==0) root_budget=GPUMemoryBudget(free_bytes,total_bytes,cap)
      END IF
      CALL MPI_Bcast(root_budget,1,MPINTLONG,0,comm,ierr)
      root_fit=root_peak<=root_budget
    END IF
    ! Avoid even creating CUDA contexts on inactive ranks when root wins.
    IF(code==2 .OR. (code==0 .AND. root_fit .AND. root_seconds<=all_seconds)) THEN
      root_only=root_fit; enabled=root_fit
    ELSE
      IF(all_index .AND. (cap<0 .OR. all_peak<=cap)) THEN
        status=device_info(device,count,identity)
        valid=0; IF(status==0) valid=1
        CALL MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,local_comm,ierr)
        CALL MPI_Comm_size(local_comm,local_size,ierr)
        CALL GET_ENVIRONMENT_VARIABLE('YUKAWA_NCCL_ALLOW_MULTINODE',multinode)
        IF(TRIM(multinode)/='1') THEN
          IF(wh%process_grid%total_processors>4 .OR. local_size/=wh%process_grid%total_processors) valid=0
        END IF
        ALLOCATE(devices(16,local_size))
        CALL MPI_Allgather(identity,16,MPI_BYTE,devices,16,MPI_BYTE,local_comm,ierr)
        DO j=1,local_size
          DO k=j+1,local_size
            IF(ALL(devices(:,j)==devices(:,k))) valid=0
          END DO
        END DO
        CALL MPI_Comm_free(local_comm,ierr); DEALLOCATE(devices)
        CALL MPI_Allreduce(valid,all_valid,1,MPI_INTEGER,MPI_MIN,comm,ierr)
        IF(all_valid==1) THEN
          status=memory_info(free_bytes,total_bytes)
          IF(status==0) all_budget=GPUMemoryBudget(free_bytes,total_bytes,cap)
          CALL MPI_Allreduce(MPI_IN_PLACE,all_budget,1,MPINTLONG,MPI_MIN,comm,ierr)
          all_fit=all_peak<=all_budget
        END IF
      END IF
      IF(all_fit) THEN
        enabled=.TRUE.
      ELSE IF(code==0 .AND. root_fit) THEN
        enabled=.TRUE.; root_only=.TRUE.
      END IF
    END IF
    IF(rank==0) THEN
      WRITE(*,'(A,4(1X,I0))') 'ntpoly_gpu_plan_peak_budget_bytes root/all:', &
        root_peak,root_budget,all_peak,all_budget
      WRITE(*,'(A,2(1X,ES14.6))') 'ntpoly_gpu_plan_estimated_seconds root/all:',root_seconds,all_seconds
      IF(.NOT.enabled) PRINT *, 'NTPoly GPU fallback: CPU (no admissible GPU layout)'
    END IF
  END SUBROUTINE

  ! Recheck after runtime/NCCL creation, before allocating matrix payloads.
  SUBROUTINE CheckGPUMemory(wh,x,root_only,allowed,extra_fixed_nnz)
    TYPE(Matrix_ps), INTENT(IN) :: wh,x
    LOGICAL, INTENT(IN) :: root_only
    LOGICAL, INTENT(OUT) :: allowed
    INTEGER(INT64), OPTIONAL, INTENT(IN) :: extra_fixed_nnz
    INTEGER :: comm,rank,ierr,j,k,valid,all_valid,pr,pc,n,io
    INTEGER(INT64) :: local_wh,nnz_wh,peak,budget,cap,cap_mb,free_bytes,total_bytes
    REAL(REAL64) :: seconds
    LOGICAL :: index_safe,participant
    CHARACTER(64) :: cap_text
    INTEGER(C_INT) :: status
    comm=wh%process_grid%global_comm; rank=wh%process_grid%global_rank
    participant=.NOT.root_only .OR. rank==0
    local_wh=0
    DO j=1,SIZE(wh%local_data_r,2)
      DO k=1,SIZE(wh%local_data_r,1)
        local_wh=local_wh+SIZE(wh%local_data_r(k,j)%values,KIND=INT64)
      END DO
    END DO
    CALL MPI_Allreduce(local_wh,nnz_wh,1,MPINTLONG,MPI_SUM,comm,ierr)
    IF(PRESENT(extra_fixed_nnz)) THEN
      IF(extra_fixed_nnz<0 .OR. nnz_wh>HUGE(nnz_wh)-extra_fixed_nnz) CALL FailGPU('fixed input count overflow')
      nnz_wh=nnz_wh+extra_fixed_nnz
    END IF
    pr=wh%process_grid%num_process_rows; pc=wh%process_grid%num_process_columns
    n=wh%logical_matrix_dimension
    IF(root_only) THEN
      pr=1; pc=1; n=wh%actual_matrix_dimension
    END IF
    CALL EstimateGPUPlan(n,pr,pc,nnz_wh,0_INT64,1,peak,seconds,index_safe)
    CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_MEMORY_MB',cap_text)
    cap=-1
    IF(LEN_TRIM(cap_text)>0) THEN
      READ(cap_text,*,IOSTAT=io) cap_mb
      IF(io/=0) CALL FailGPU('invalid GPU memory cap')
      cap=cap_mb*1048576_INT64
    END IF
    valid=1
    IF(participant) THEN
      status=memory_info(free_bytes,total_bytes)
      budget=0
      IF(status==0) budget=GPUMemoryBudget(free_bytes,total_bytes,cap)
      IF(.NOT.index_safe .OR. peak>budget) valid=0
    END IF
    CALL MPI_Allreduce(valid,all_valid,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    allowed=all_valid==1
    IF(.NOT.allowed .AND. rank==0) PRINT *, 'NTPoly GPU fallback: CPU (post-bootstrap memory admission)'
  END SUBROUTINE

  SUBROUTINE BeginGPUBackend(matrix,root_only)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    LOGICAL, OPTIONAL, INTENT(IN) :: root_only
    INTEGER(C_INT) :: status,device
    INTEGER :: mode,comm,hit,all_hit,ierr
    LOGICAL :: root_mode,reused
    root_mode=.FALSE.
    IF(PRESENT(root_only)) root_mode=root_only
    IF(root_mode .AND. matrix%process_grid%global_rank/=0) RETURN
    IF(nccl_active) CALL FailGPU('nested GPU solver contexts are unsupported')
    CALL init_gpu(status,device)
    IF(status/=0) CALL FailGPU('CUDA initialization failed')
    mode=1; comm=matrix%process_grid%global_comm
    IF(root_mode) THEN
      mode=2; comm=MPI_COMM_SELF
    END IF
    reused=C_ASSOCIATED(matrix%process_grid%gpu_context(mode))
    hit=0; IF(reused) hit=1
    CALL MPI_Allreduce(hit,all_hit,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    CALL MPI_Comm_size(comm,hit,ierr)
    IF(all_hit/=0 .AND. all_hit/=hit) CALL FailGPU('inconsistent grid context ownership')
    IF(.NOT.reused) THEN
      IF(create_context(matrix%process_grid%gpu_context(mode))/=0) CALL FailGPU('context allocation failed')
    END IF
    IF(select_context(matrix%process_grid%gpu_context(mode))/=0) CALL FailGPU('context activation failed')
    IF(reused) THEN
      nccl_active=.TRUE.
    ELSE
      IF(root_mode) THEN
        CALL InitSolverNCCL(MPI_COMM_SELF,MPI_COMM_SELF,MPI_COMM_SELF)
        CALL WarmSolverNCCL(MPI_COMM_SELF,MPI_COMM_SELF)
      ELSE
        CALL InitSolverNCCL(matrix%process_grid%row_comm,matrix%process_grid%column_comm,matrix%process_grid%global_comm)
        CALL WarmSolverNCCL(matrix%process_grid%row_comm,matrix%process_grid%column_comm)
      END IF
    END IF
    IF(matrix%process_grid%global_rank==0) THEN
      PRINT *, 'NTPoly GPU context reused:',reused
      IF(root_mode) THEN
        PRINT *, 'NTPoly GPU TRS2 active ranks: 1 (coalesced)'
      ELSE
        PRINT *, 'NTPoly GPU TRS2 active ranks:',matrix%process_grid%total_processors
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE EndGPUBackend()
    IF(free_panels()/=0) CALL FailGPU('panel cleanup failed')
    IF(suspend_context()/=0) CALL FailGPU('context release failed')
    nccl_active=.FALSE.
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE DestroyGPU(matrix)
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    IF(matrix%handle/=0) CALL destroy_slot(matrix%handle)
  END SUBROUTINE

  SUBROUTINE ReleaseGPU(matrix)
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    CALL DestroyGPU(matrix)
  END SUBROUTINE

  SUBROUTINE ShapeGPU(host,matrix)
    TYPE(Matrix_ps), INTENT(IN) :: host
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    IF(matrix%handle==0) CALL create_slot(matrix%handle)
    matrix%rows=host%local_rows; matrix%columns=host%local_columns
    matrix%dimension=host%logical_matrix_dimension
    matrix%diagonal=host%start_row-host%start_column
    matrix%process_rows=host%process_grid%num_process_rows
    matrix%process_columns=host%process_grid%num_process_columns
  END SUBROUTINE

  SUBROUTINE UploadGPU(host,matrix)
    TYPE(Matrix_ps), INTENT(IN), TARGET :: host
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    TYPE(Matrix_lsr) :: local
    IF(host%is_complex) CALL FailGPU('complex upload is unsupported')
    CALL ShapeGPU(host,matrix)
    IF(SIZE(host%local_data_r,1)==1 .AND. SIZE(host%local_data_r,2)==1) THEN
      ! Borrow the existing block for the synchronous conversion. No merge
      ! or payload access on the CPU is needed for the common one-block case.
      CALL UploadLocalGPU(host%local_data_r(1,1),matrix)
      RETURN
    END IF
    CALL MergeMatrixLocalBlocks(host,local)
    CALL UploadLocalGPU(local,matrix)
    CALL DestructMatrix(local)
  END SUBROUTINE

  SUBROUTINE UploadLocalGPU(local,matrix)
    TYPE(Matrix_lsr), INTENT(IN), TARGET :: local
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    INTEGER, ALLOCATABLE :: rp(:),ci(:)
    REAL(C_DOUBLE), ALLOCATABLE :: values(:)
    INTEGER(C_INT) :: nnz,status
    TYPE(C_PTR) :: pr,pi,pv
    INTEGER(C_INT), POINTER :: fr(:),fi(:)
    REAL(C_DOUBLE), POINTER :: fv(:)
    CHARACTER(16) :: mode,report
    mode=''
    CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_IMPORT',mode)
    IF(TRIM(mode)/='' .AND. TRIM(mode)/='auto' .AND. TRIM(mode)/='host') CALL FailGPU('invalid GPU import mode')
    IF(TRIM(mode)/='host') THEN
      nnz=SIZE(local%values)
      pi=C_NULL_PTR; pv=C_NULL_PTR
      IF(nnz>0) THEN
        pi=C_LOC(local%inner_index(1)); pv=C_LOC(local%values(1))
      END IF
      status=import_csc(matrix%handle,matrix%rows,matrix%columns,nnz,C_LOC(local%outer_index(1)),pi,pv)
      IF(status<0) CALL FailGPU('device CSC import failed')
      IF(status==0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_IMPORT_REPORT',report)
        IF(TRIM(report)=='1') PRINT *, 'NTPoly GPU import: borrowed CSC arrays, nnz=',nnz
        RETURN
      END IF
    END IF
    CALL TransposeCSR(matrix%columns,matrix%rows,local%outer_index,local%inner_index,local%values,rp,ci,values)
    nnz=SIZE(ci)
    CALL allocate_slot(matrix%handle,matrix%rows,nnz,pr,pi,pv,status)
    IF(status/=0) CALL FailGPU('matrix allocation failed')
    CALL C_F_POINTER(pr,fr,[matrix%rows+1]); CALL C_F_POINTER(pi,fi,[MAX(nnz,1)])
    CALL C_F_POINTER(pv,fv,[MAX(nnz,1)])
    fr=rp; fi(1:nnz)=ci; fv(1:nnz)=values
  END SUBROUTINE

  SUBROUTINE DownloadLocalGPU(matrix,local)
    TYPE(GPUMatrix_t), INTENT(IN) :: matrix
    TYPE(Matrix_lsr), INTENT(INOUT) :: local
    INTEGER(C_INT) :: rows,nnz
    TYPE(C_PTR) :: pr,pi,pv
    INTEGER(C_INT), POINTER :: rp(:),ci(:)
    REAL(C_DOUBLE), POINTER :: values(:)
    IF(matrix%handle==0) CALL FailGPU('download from uninitialized matrix')
    CALL get_slot(matrix%handle,pr,pi,pv,rows,nnz)
    CALL C_F_POINTER(pr,rp,[rows+1]); CALL C_F_POINTER(pi,ci,[MAX(nnz,1)])
    CALL C_F_POINTER(pv,values,[MAX(nnz,1)])
    local%rows=matrix%rows; local%columns=matrix%columns
    CALL TransposeCSR(matrix%rows,matrix%columns,rp,ci,values,local%outer_index,local%inner_index,local%values)
  END SUBROUTINE

  SUBROUTINE DownloadGPU(matrix,host)
    TYPE(GPUMatrix_t), INTENT(IN) :: matrix
    TYPE(Matrix_ps), INTENT(INOUT) :: host
    TYPE(Matrix_lsr) :: local
    IF(host%is_complex .OR. host%local_rows/=matrix%rows .OR. host%local_columns/=matrix%columns) &
      CALL FailGPU('download shape mismatch')
    CALL DownloadLocalGPU(matrix,local)
    CALL SplitMatrixToLocalBlocks(host,local)
    CALL DestructMatrix(local)
  END SUBROUTINE

  SUBROUTINE UploadGPURegion(host,matrix,root_only)
    TYPE(Matrix_ps), INTENT(IN) :: host
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    LOGICAL, INTENT(IN) :: root_only
    TYPE(Matrix_lsr) :: local
    IF(.NOT.root_only) THEN
      CALL UploadGPU(host,matrix)
      RETURN
    END IF
    CALL GatherMatrixToProcess(host,local,0)
    IF(host%process_grid%global_rank==0) THEN
      IF(matrix%handle==0) CALL create_slot(matrix%handle)
      matrix%rows=host%actual_matrix_dimension; matrix%columns=matrix%rows
      matrix%dimension=matrix%rows; matrix%diagonal=0
      matrix%process_rows=1; matrix%process_columns=1
      CALL UploadLocalGPU(local,matrix)
    END IF
    CALL DestructMatrix(local)
  END SUBROUTINE

  SUBROUTINE DownloadGPURegion(matrix,host,root_only)
    TYPE(GPUMatrix_t), INTENT(IN) :: matrix
    TYPE(Matrix_ps), INTENT(INOUT) :: host
    LOGICAL, INTENT(IN) :: root_only
    TYPE(Matrix_lsr) :: local
    TYPE(TripletList_r) :: entries
    IF(.NOT.root_only) THEN
      CALL DownloadGPU(matrix,host)
      RETURN
    END IF
    CALL ConstructTripletList(entries)
    IF(host%process_grid%global_rank==0) THEN
      CALL DownloadLocalGPU(matrix,local)
      CALL MatrixToTripletList(local,entries)
    END IF
    CALL FillMatrixFromTripletList(host,entries,preduplicated_in=.TRUE.)
    CALL DestructMatrix(local); CALL DestructTripletList(entries)
  END SUBROUTINE

  SUBROUTINE CopyGPU(source,destination)
    CLASS(GPUMatrix_t), INTENT(IN) :: source
    CLASS(GPUMatrix_t), INTENT(INOUT) :: destination
    IF(source%handle==destination%handle .AND. source%handle/=0) RETURN
    IF(source%handle==0) CALL FailGPU('copy from uninitialized matrix')
    IF(destination%handle==0) CALL create_slot(destination%handle)
    CALL copy_slot(destination%handle,source%handle)
    destination%rows=source%rows; destination%columns=source%columns
    destination%dimension=source%dimension; destination%diagonal=source%diagonal
    destination%process_rows=source%process_rows; destination%process_columns=source%process_columns
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE AssignGPU(destination,source)
    CLASS(GPUMatrix_t), INTENT(INOUT) :: destination
    TYPE(GPUMatrix_t), INTENT(IN) :: source
    IF(source%handle==0) THEN
      IF(destination%handle/=0) CALL destroy_slot(destination%handle)
    ELSE
      CALL CopyGPU(source,destination)
    END IF
  END SUBROUTINE

  SUBROUTINE ScaleGPU(matrix,alpha)
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    REAL(NTREAL), INTENT(IN) :: alpha
    IF(matrix%handle==0) CALL FailGPU('scale of uninitialized matrix')
    CALL scale_slot(matrix%handle,alpha)
  END SUBROUTINE

  SUBROUTINE CheckShapes(a,b)
    TYPE(GPUMatrix_t), INTENT(IN) :: a,b
    IF(a%handle==0 .OR. b%handle==0) CALL FailGPU('operation on uninitialized matrix')
    IF(a%rows/=b%rows .OR. a%columns/=b%columns .OR. a%dimension/=b%dimension .OR. &
       a%process_rows/=b%process_rows .OR. a%process_columns/=b%process_columns .OR. a%diagonal/=b%diagonal) &
      CALL FailGPU('matrix layouts differ')
  END SUBROUTINE

  SUBROUTINE IncrementGPU(source,destination,alpha,threshold)
    TYPE(GPUMatrix_t), INTENT(IN) :: source
    TYPE(GPUMatrix_t), INTENT(INOUT) :: destination
    REAL(NTREAL), INTENT(IN) :: alpha,threshold
    INTEGER(C_INT) :: nnz
    CALL CheckShapes(source,destination)
    CALL increment_slot(destination%handle,source%handle,alpha,threshold,nnz)
    IF(nnz<0) CALL FailGPU('increment failed')
  END SUBROUTINE

  SUBROUTINE TraceGPU(matrix,value)
    TYPE(GPUMatrix_t), INTENT(IN) :: matrix
    REAL(NTREAL), INTENT(OUT) :: value
    IF(matrix%handle==0) CALL FailGPU('trace of uninitialized matrix')
    CALL trace_slot(matrix%handle,matrix%diagonal,value)
    CALL SolverScalarSum(value)
  END SUBROUTINE

  SUBROUTINE BoundsGPURegion(matrix,host,root_only,lower_bound,upper_bound)
    TYPE(GPUMatrix_t), INTENT(IN) :: matrix
    TYPE(Matrix_ps), INTENT(IN) :: host
    LOGICAL, INTENT(IN) :: root_only
    REAL(NTREAL), INTENT(OUT) :: lower_bound,upper_bound
    REAL(C_DOUBLE), ALLOCATABLE :: lower(:),upper(:)
    REAL(NTREAL) :: bounds(2)
    INTEGER :: ierr
    IF(.NOT.root_only .OR. host%process_grid%global_rank==0) THEN
      ALLOCATE(lower(matrix%columns),upper(matrix%columns))
      IF(bounds_slot(matrix%handle,matrix%columns,matrix%diagonal,lower,upper)/=0) CALL FailGPU('bounds failed')
      IF(.NOT.root_only) THEN
        CALL MPI_Allreduce(MPI_IN_PLACE,lower,SIZE(lower),MPINTREAL,MPI_SUM,host%process_grid%column_comm,ierr)
        CALL MPI_Allreduce(MPI_IN_PLACE,upper,SIZE(upper),MPINTREAL,MPI_SUM,host%process_grid%column_comm,ierr)
      END IF
      lower_bound=MINVAL(lower); upper_bound=MAXVAL(upper)
      IF(root_only) THEN
        ! The coalesced CSR omits padded zero columns; the CPU algorithm
        ! includes them in its extrema. Retain that numerical convention.
        IF(host%logical_matrix_dimension>host%actual_matrix_dimension) THEN
          lower_bound=MIN(0.0_NTREAL,lower_bound); upper_bound=MAX(0.0_NTREAL,upper_bound)
        END IF
      ELSE
        CALL MPI_Allreduce(MPI_IN_PLACE,lower_bound,1,MPINTREAL,MPI_MIN,host%process_grid%row_comm,ierr)
        CALL MPI_Allreduce(MPI_IN_PLACE,upper_bound,1,MPINTREAL,MPI_MAX,host%process_grid%row_comm,ierr)
      END IF
      DEALLOCATE(lower,upper)
      bounds=[lower_bound,upper_bound]
    END IF
    IF(root_only) THEN
      CALL MPI_Bcast(bounds,2,MPINTREAL,0,host%process_grid%global_comm,ierr)
      lower_bound=bounds(1); upper_bound=bounds(2)
    END IF
  END SUBROUTINE

  SUBROUTINE DotGPU(a,b,value)
    TYPE(GPUMatrix_t), INTENT(IN) :: a,b
    REAL(NTREAL), INTENT(OUT) :: value
    CALL CheckShapes(a,b)
    CALL dot_slot(a%handle,b%handle,value)
    CALL SolverScalarSum(value)
  END SUBROUTINE

  SUBROUTINE SquareGPU(source,destination,threshold)
    TYPE(GPUMatrix_t), INTENT(IN) :: source
    TYPE(GPUMatrix_t), INTENT(INOUT) :: destination
    REAL(NTREAL), INTENT(IN) :: threshold
    CALL MultiplyGPU(source,source,destination,threshold)
  END SUBROUTINE

  SUBROUTINE MultiplyGPU(left,right,destination,threshold)
    TYPE(GPUMatrix_t), INTENT(IN) :: left,right
    TYPE(GPUMatrix_t), INTENT(INOUT) :: destination
    REAL(NTREAL), INTENT(IN) :: threshold
    REAL(C_DOUBLE) :: phases(4)
    CALL CheckShapes(left,right)
    IF(destination%handle==0) CALL create_slot(destination%handle)
    IF(left%handle==destination%handle .OR. right%handle==destination%handle) CALL FailGPU('in-place product is unsupported')
    IF(resident_multiply(left%handle,right%handle,destination%handle,left%rows,left%dimension,left%columns, &
        left%process_rows,left%process_columns,threshold,phases)/=0) CALL FailGPU('resident product failed')
    destination%rows=left%rows; destination%columns=left%columns
    destination%dimension=left%dimension; destination%diagonal=left%diagonal
    destination%process_rows=left%process_rows; destination%process_columns=left%process_columns
  END SUBROUTINE

  ! Transitional host-facing transform: products are resident, boundaries still
  ! import/export Matrix_ps. Preserve the CPU path when admission is unsafe.
  SUBROUTINE TryGPUSimilarityTransform(a,p,pinv,result,threshold,performed)
    TYPE(Matrix_ps), INTENT(IN) :: a,p,pinv
    TYPE(Matrix_ps), INTENT(INOUT) :: result
    REAL(NTREAL), INTENT(IN) :: threshold
    LOGICAL, INTENT(OUT) :: performed
    TYPE(GPUMatrix_t) :: ga,gp,gpinv,temp,output
    LOGICAL :: root_mode,allowed
    INTEGER :: j,k,ierr,valid,all_valid
    INTEGER(INT64) :: extra,global_extra
    performed=.FALSE.; valid=1
    IF(a%is_complex .OR. p%is_complex .OR. pinv%is_complex) valid=0
    IF(.NOT.ASSOCIATED(a%process_grid,p%process_grid) .OR. &
       .NOT.ASSOCIATED(a%process_grid,pinv%process_grid)) valid=0
    IF(a%logical_matrix_dimension/=p%logical_matrix_dimension .OR. &
       a%logical_matrix_dimension/=pinv%logical_matrix_dimension) valid=0
    IF(a%actual_matrix_dimension/=p%actual_matrix_dimension .OR. &
       a%actual_matrix_dimension/=pinv%actual_matrix_dimension) valid=0
    CALL MPI_Allreduce(valid,all_valid,1,MPI_INTEGER,MPI_MIN,a%process_grid%global_comm,ierr)
    IF(all_valid==0) RETURN
    ! Three fixed inputs, two product outputs; the planner's three full dynamic
    ! CSR slots conservatively cover both outputs plus its existing scratch.
    extra=0
    DO j=1,SIZE(p%local_data_r,2)
      DO k=1,SIZE(p%local_data_r,1)
        extra=extra+SIZE(p%local_data_r(k,j)%values,KIND=INT64)
        extra=extra+SIZE(pinv%local_data_r(k,j)%values,KIND=INT64)
      END DO
    END DO
    CALL MPI_Allreduce(extra,global_extra,1,MPINTLONG,MPI_SUM,a%process_grid%global_comm,ierr)
    CALL SelectGPUBackend(a,a,allowed,root_mode,2,global_extra)
    IF(.NOT.allowed) RETURN
    CALL BeginGPUBackend(a,root_mode)
    CALL CheckGPUMemory(a,a,root_mode,allowed,global_extra)
    IF(allowed) THEN
      CALL UploadGPURegion(a,ga,root_mode)
      CALL UploadGPURegion(p,gp,root_mode)
      CALL UploadGPURegion(pinv,gpinv,root_mode)
      CALL ConstructEmptyMatrix(result,a)
      IF(.NOT.root_mode .OR. a%process_grid%global_rank==0) THEN
        CALL MultiplyGPU(gp,ga,temp,threshold)
        CALL MultiplyGPU(temp,gpinv,output,threshold)
      END IF
      CALL DownloadGPURegion(output,result,root_mode)
      performed=.TRUE.
    END IF
    IF(.NOT.root_mode .OR. a%process_grid%global_rank==0) CALL EndGPUBackend()
  END SUBROUTINE
END MODULE
