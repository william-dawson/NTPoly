! Experimental owned GPU matrices for the library TRS2 iteration.
! Host Matrix_ps remains authoritative outside the explicitly bounded GPU region.
MODULE PSMatrixGPUBackendModule
  USE ISO_C_BINDING
  USE DataTypesModule, ONLY: NTREAL, MPINTLONG, MPINTREAL
  USE ISO_FORTRAN_ENV, ONLY: INT64, REAL64
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE GPUPlanningModule, ONLY: EstimateGPUPlan, EstimateComplexProductPlan, GPUMemoryBudget
  USE NTMPIModule
  USE PSMatrixModule, ONLY: Matrix_ps, ConstructEmptyMatrix, MergeMatrixLocalBlocks, SplitMatrixToLocalBlocks, &
    GatherMatrixToProcess, FillMatrixFromTripletList, GetMatrixTripletList
  USE SMatrixModule, ONLY: Matrix_lsr, Matrix_lsc, DestructMatrix, MatrixToTripletList
  USE TripletListModule, ONLY: TripletList_r, TripletList_c, ConstructTripletList, DestructTripletList
  USE SolverCSR, ONLY: TransposeCSR
  USE SolverNCCL, ONLY: InitSolverNCCL, WarmSolverNCCL, FinishSolverNCCL, SolverScalarSum, nccl_active
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GPUMatrix_t, GPURequested, SelectGPUBackend, CheckGPUMemory, BeginGPUBackend, EndGPUBackend
  PUBLIC :: UploadGPURegion, DownloadGPURegion
  PUBLIC :: UploadGPU, DownloadGPU, CopyGPU, ScaleGPU, IncrementGPU, TraceGPU, DotGPU, SquareGPU, MultiplyGPU
  PUBLIC :: TryGPUSimilarityTransform, TryGPUMatrixMultiply, ReleaseGPU, BoundsGPURegion, PruneGPU
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
    INTEGER(C_INT) FUNCTION complex_product(rows,columns,n,row_offset,column_offset,na,ar,ac,av,nb,br,bc,bv, &
        nc,cr,cc,cv,alpha,beta,threshold,output,count) BIND(C,NAME='ntpoly_complex_product_')
      IMPORT C_INT,C_DOUBLE,C_DOUBLE_COMPLEX,C_PTR
      INTEGER(C_INT) :: rows,columns,n,row_offset,column_offset,na,nb,nc,count
      INTEGER(C_INT) :: ar(*),ac(*),br(*),bc(*),cr(*),cc(*)
      COMPLEX(C_DOUBLE_COMPLEX) :: av(*),bv(*),cv(*)
      REAL(C_DOUBLE) :: alpha,beta,threshold
      TYPE(C_PTR) :: output
    END FUNCTION
    SUBROUTINE complex_export(output,rows,columns,values) BIND(C,NAME='ntpoly_complex_product_export_')
      IMPORT C_INT,C_DOUBLE_COMPLEX,C_PTR
      TYPE(C_PTR) :: output
      INTEGER(C_INT) :: rows(*),columns(*)
      COMPLEX(C_DOUBLE_COMPLEX) :: values(*)
    END SUBROUTINE
    INTEGER(C_INT) FUNCTION clear_context_workspace(context,limit,cleared) &
        BIND(C,NAME='solver_nccl_context_clear_workspace_')
      IMPORT C_PTR,C_INT,C_INT64_T
      TYPE(C_PTR), VALUE :: context
      INTEGER(C_INT64_T) :: limit
      INTEGER(C_INT) :: cleared
    END FUNCTION
    INTEGER(C_INT) FUNCTION clear_active_workspace() BIND(C,NAME='solver_nccl_workspace_clear_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT) FUNCTION trim_workspaces(cleared) BIND(C,NAME='solver_nccl_trim_workspaces_')
      IMPORT C_INT
      INTEGER(C_INT) :: cleared
    END FUNCTION
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
    SUBROUTINE prune_slot(slot,threshold,nnz) BIND(C,NAME='gpu_ops_prune_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT) :: slot,nnz
      REAL(C_DOUBLE) :: threshold
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

  LOGICAL FUNCTION GPURequested(matrix,option_name) RESULT(requested)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    CHARACTER(*), OPTIONAL, INTENT(IN) :: option_name
    CHARACTER(32) :: name
    CHARACTER(16) :: option
    INTEGER :: flags(2),global_flags(2),ierr
    name='NTPOLY_GPU_TRS2'
    IF(PRESENT(option_name)) name=option_name
    CALL GET_ENVIRONMENT_VARIABLE(TRIM(name),option)
    flags(1)=0
    IF(TRIM(option)=='1') flags(1)=1
    flags(2)=-flags(1)
    CALL MPI_Allreduce(flags,global_flags,2,MPI_INTEGER,MPI_MAX,matrix%process_grid%global_comm,ierr)
    IF(global_flags(1)/=-global_flags(2)) CALL FailGPU(TRIM(name)//' differs between ranks')
    requested=global_flags(1)==1
    IF(.NOT.requested) RETURN
    IF(matrix%is_complex .AND. .NOT.PRESENT(option_name)) THEN
      ! Complex TRS2 can use shared products, but the resident real iteration
      ! must not reinterpret complex payloads or abort an otherwise valid solve.
      requested=.FALSE.;RETURN
    END IF
    IF(matrix%process_grid%num_process_slices/=1) THEN
      IF(PRESENT(option_name)) THEN
        requested=.FALSE.
      ELSE
        CALL FailGPU('TRS2 GPU backend currently requires real matrices and one slice')
      END IF
    END IF
  END FUNCTION

  ! Selection happens after CPU WH/X0 setup, using the actual working pattern.
  RECURSIVE SUBROUTINE SelectGPUBackend(wh,x,enabled,root_only,expected_iterations,extra_fixed_nnz, &
      initial_shift,preserve_padding,memory_limit,reclaim_workspace)
    TYPE(Matrix_ps), INTENT(IN) :: wh,x
    LOGICAL, INTENT(OUT) :: enabled,root_only
    INTEGER, INTENT(IN) :: expected_iterations
    INTEGER(INT64), OPTIONAL, INTENT(IN) :: extra_fixed_nnz
    LOGICAL, OPTIONAL, INTENT(IN) :: initial_shift,preserve_padding
    INTEGER(INT64), OPTIONAL, INTENT(OUT) :: memory_limit
    LOGICAL, OPTIONAL, INTENT(IN) :: reclaim_workspace
    CHARACTER(64) :: policy,cap_text,multinode
    INTEGER :: code,flags(2),global_flags(2),ierr,pr,pc,comm,rank,iterations,io,j,k
    INTEGER :: local_comm,local_size,valid,all_valid,root_dimension
    INTEGER(C_INT8_T), ALLOCATABLE :: devices(:,:)
    INTEGER(C_INT8_T) :: identity(16)
    INTEGER(C_INT) :: device,count,status
    INTEGER(INT64) :: nnz_wh,nnz_x,root_peak,all_peak,root_budget,all_budget,cap,local_cap
    INTEGER(INT64) :: cap_mb,local_wh,local_x,free_bytes,total_bytes
    REAL(REAL64) :: root_seconds,all_seconds
    LOGICAL :: root_index,all_index,root_fit,all_fit,reclaim
    comm=wh%process_grid%global_comm; rank=wh%process_grid%global_rank
    enabled=.FALSE.; root_only=.FALSE.
    reclaim=.TRUE.;IF(PRESENT(reclaim_workspace)) reclaim=reclaim_workspace
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
    IF(PRESENT(memory_limit)) memory_limit=cap
    ! Count local sub-blocks without merging/copying their payloads.
    local_wh=0; local_x=0
    local_wh=LocalGPUCount(wh);local_x=LocalGPUCount(x)
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
    root_dimension=wh%actual_matrix_dimension
    IF(PRESENT(preserve_padding)) THEN
      IF(preserve_padding) root_dimension=wh%logical_matrix_dimension
    END IF
    CALL MatrixGPUPlan(wh,root_dimension,1,1,nnz_wh,nnz_x,iterations,root_peak,root_seconds,root_index)
    CALL MatrixGPUPlan(wh,wh%logical_matrix_dimension,pr,pc,nnz_wh,nnz_x,iterations,all_peak,all_seconds,all_index)
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
    IF(.NOT.enabled .AND. reclaim) THEN
      IF(EvictCachedGPUWorkspaces(wh)) THEN
        CALL SelectGPUBackend(wh,x,enabled,root_only,expected_iterations,extra_fixed_nnz, &
          initial_shift,preserve_padding,memory_limit,reclaim_workspace=.FALSE.)
        RETURN
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
  SUBROUTINE CheckGPUMemory(wh,x,root_only,allowed,extra_fixed_nnz,preserve_padding)
    TYPE(Matrix_ps), INTENT(IN) :: wh,x
    LOGICAL, INTENT(IN) :: root_only
    LOGICAL, INTENT(OUT) :: allowed
    INTEGER(INT64), OPTIONAL, INTENT(IN) :: extra_fixed_nnz
    LOGICAL, OPTIONAL, INTENT(IN) :: preserve_padding
    INTEGER :: comm,rank,ierr,j,k,valid,all_valid,pr,pc,n,io
    INTEGER(INT64) :: local_wh,nnz_wh,peak,budget,cap,cap_mb,free_bytes,total_bytes
    REAL(REAL64) :: seconds
    LOGICAL :: index_safe,participant
    CHARACTER(64) :: cap_text
    INTEGER(C_INT) :: status
    comm=wh%process_grid%global_comm; rank=wh%process_grid%global_rank
    participant=.NOT.root_only .OR. rank==0
    local_wh=0
    local_wh=LocalGPUCount(wh)
    CALL MPI_Allreduce(local_wh,nnz_wh,1,MPINTLONG,MPI_SUM,comm,ierr)
    IF(PRESENT(extra_fixed_nnz)) THEN
      IF(extra_fixed_nnz<0 .OR. nnz_wh>HUGE(nnz_wh)-extra_fixed_nnz) CALL FailGPU('fixed input count overflow')
      nnz_wh=nnz_wh+extra_fixed_nnz
    END IF
    pr=wh%process_grid%num_process_rows; pc=wh%process_grid%num_process_columns
    n=wh%logical_matrix_dimension
    IF(root_only) THEN
      pr=1; pc=1; n=wh%actual_matrix_dimension
      IF(PRESENT(preserve_padding)) THEN
        IF(preserve_padding) n=wh%logical_matrix_dimension
      END IF
    END IF
    CALL MatrixGPUPlan(wh,n,pr,pc,nnz_wh,0_INT64,1,peak,seconds,index_safe)
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

  SUBROUTINE UploadGPURegion(host,matrix,root_only,preserve_padding)
    TYPE(Matrix_ps), INTENT(IN) :: host
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    LOGICAL, INTENT(IN) :: root_only
    LOGICAL, OPTIONAL, INTENT(IN) :: preserve_padding
    LOGICAL :: keep_padding
    TYPE(Matrix_lsr) :: local
    IF(.NOT.root_only) THEN
      CALL UploadGPU(host,matrix)
      RETURN
    END IF
    keep_padding=.FALSE.
    IF(PRESENT(preserve_padding)) keep_padding=preserve_padding
    CALL GatherMatrixToProcess(host,local,0,include_padding_in=keep_padding)
    IF(host%process_grid%global_rank==0) THEN
      IF(matrix%handle==0) CALL create_slot(matrix%handle)
      matrix%rows=host%actual_matrix_dimension
      IF(keep_padding) matrix%rows=host%logical_matrix_dimension
      matrix%columns=matrix%rows
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

  SUBROUTINE PruneGPU(matrix,threshold)
    TYPE(GPUMatrix_t), INTENT(INOUT) :: matrix
    REAL(NTREAL), INTENT(IN) :: threshold
    INTEGER(C_INT) :: nnz
    IF(matrix%handle==0) CALL FailGPU('prune of uninitialized matrix')
    CALL prune_slot(matrix%handle,threshold,nnz)
    IF(nnz<0) CALL FailGPU('prune failed')
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

  SUBROUTINE MultiplyGPU(left,right,destination,threshold,phase_times)
    TYPE(GPUMatrix_t), INTENT(IN) :: left,right
    TYPE(GPUMatrix_t), INTENT(INOUT) :: destination
    REAL(NTREAL), INTENT(IN) :: threshold
    REAL(C_DOUBLE), OPTIONAL, INTENT(OUT) :: phase_times(4)
    REAL(C_DOUBLE) :: phases(4)
    CALL CheckShapes(left,right)
    IF(destination%handle==0) CALL create_slot(destination%handle)
    IF(left%handle==destination%handle .OR. right%handle==destination%handle) CALL FailGPU('in-place product is unsupported')
    IF(resident_multiply(left%handle,right%handle,destination%handle,left%rows,left%dimension,left%columns, &
        left%process_rows,left%process_columns,threshold,phases)/=0) CALL FailGPU('resident product failed')
    IF(PRESENT(phase_times)) phase_times=phases
    destination%rows=left%rows; destination%columns=left%columns
    destination%dimension=left%dimension; destination%diagonal=left%diagonal
    destination%process_rows=left%process_rows; destination%process_columns=left%process_columns
  END SUBROUTINE

  ! Shared API boundary for iterative solver families. Keep the public CPU
  ! representation and perform the expensive product on the admitted GPU grid.
  SUBROUTINE TryGPUMatrixMultiply(a,b,c,alpha,beta,threshold,performed)
    TYPE(Matrix_ps), INTENT(IN) :: a,b
    TYPE(Matrix_ps), INTENT(INOUT) :: c
    REAL(NTREAL), INTENT(IN) :: alpha,beta,threshold
    LOGICAL, INTENT(OUT) :: performed
    TYPE(GPUMatrix_t) :: ga,gb,gc,product
    LOGICAL :: root_mode,allowed,keep_c,participant
    INTEGER :: valid,all_valid,ierr,j,k
    INTEGER(INT64) :: extra,global_extra
    CHARACTER(16) :: report
    performed=.FALSE.
    IF(.NOT.GPURequested(a,'NTPOLY_GPU_ALGEBRA')) RETURN
    IF(a%is_complex .OR. b%is_complex) THEN
      CALL TryComplexMultiply(a,b,c,alpha,beta,threshold,performed)
      RETURN
    END IF
    keep_c=ABS(beta)>=TINY(beta)
    valid=1
    IF(nccl_active .OR. b%is_complex) valid=0
    IF(.NOT.ASSOCIATED(a%process_grid,b%process_grid)) valid=0
    IF(a%logical_matrix_dimension/=b%logical_matrix_dimension .OR. &
       a%actual_matrix_dimension/=b%actual_matrix_dimension) valid=0
    IF(.NOT.IEEE_IS_FINITE(alpha) .OR. .NOT.IEEE_IS_FINITE(beta) .OR. &
       .NOT.IEEE_IS_FINITE(threshold) .OR. threshold<0) valid=0
    IF(keep_c) THEN
      IF(c%is_complex .OR. .NOT.ASSOCIATED(a%process_grid,c%process_grid)) valid=0
      IF(a%logical_matrix_dimension/=c%logical_matrix_dimension .OR. &
         a%actual_matrix_dimension/=c%actual_matrix_dimension) valid=0
    END IF
    CALL MPI_Allreduce(valid,all_valid,1,MPI_INTEGER,MPI_MIN,a%process_grid%global_comm,ierr)
    IF(all_valid==0) RETURN
    extra=0
    DO j=1,SIZE(b%local_data_r,2)
      DO k=1,SIZE(b%local_data_r,1)
        extra=extra+SIZE(b%local_data_r(k,j)%values,KIND=INT64)
      END DO
    END DO
    IF(keep_c) THEN
      DO j=1,SIZE(c%local_data_r,2)
        DO k=1,SIZE(c%local_data_r,1)
          extra=extra+SIZE(c%local_data_r(k,j)%values,KIND=INT64)
        END DO
      END DO
    END IF
    CALL MPI_Allreduce(extra,global_extra,1,MPINTLONG,MPI_SUM,a%process_grid%global_comm,ierr)
    CALL SelectGPUBackend(a,a,allowed,root_mode,1,global_extra,preserve_padding=.TRUE.)
    IF(.NOT.allowed) RETURN
    CALL BeginGPUBackend(a,root_mode)
    CALL CheckGPUMemory(a,a,root_mode,allowed,global_extra,preserve_padding=.TRUE.)
    participant=.NOT.root_mode .OR. a%process_grid%global_rank==0
    IF(allowed) THEN
      CALL UploadGPURegion(a,ga,root_mode,preserve_padding=.TRUE.)
      CALL UploadGPURegion(b,gb,root_mode,preserve_padding=.TRUE.)
      IF(keep_c) CALL UploadGPURegion(c,gc,root_mode,preserve_padding=.TRUE.)
      IF(participant) THEN
        ! CPU GEMM thresholds alpha*A*B BEFORE adding beta*C.
        ! Thresholding AB first would lose terms when |alpha|>1.
        IF(ABS(alpha)==1.0_NTREAL) THEN
          CALL MultiplyGPU(ga,gb,product,threshold)
          IF(alpha<0) CALL ScaleGPU(product,alpha)
        ELSE
          CALL MultiplyGPU(ga,gb,product,0.0_NTREAL)
          CALL ScaleGPU(product,alpha)
          CALL PruneGPU(product,threshold)
        END IF
        IF(keep_c) THEN
          CALL ScaleGPU(gc,beta)
          CALL IncrementGPU(product,gc,1.0_NTREAL,0.0_NTREAL)
        END IF
      END IF
      CALL ConstructEmptyMatrix(c,a)
      IF(keep_c) THEN
        CALL DownloadGPURegion(gc,c,root_mode)
      ELSE
        CALL DownloadGPURegion(product,c,root_mode)
      END IF
      performed=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_ALGEBRA_REPORT',report)
      IF(TRIM(report)=='1' .AND. a%process_grid%global_rank==0) PRINT *, 'NTPoly GPU algebra product'
    END IF
    CALL ReleaseGPU(ga); CALL ReleaseGPU(gb); CALL ReleaseGPU(gc); CALL ReleaseGPU(product)
    IF(participant) CALL EndGPUBackend()
  END SUBROUTINE

  INTEGER(INT64) FUNCTION LocalGPUCount(matrix) RESULT(count)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    INTEGER :: i,j
    count=0
    IF(matrix%is_complex) THEN
      DO j=1,SIZE(matrix%local_data_c,2)
        DO i=1,SIZE(matrix%local_data_c,1)
          count=count+SIZE(matrix%local_data_c(i,j)%values,KIND=INT64)
        END DO
      END DO
    ELSE
      DO j=1,SIZE(matrix%local_data_r,2)
        DO i=1,SIZE(matrix%local_data_r,1)
          count=count+SIZE(matrix%local_data_r(i,j)%values,KIND=INT64)
        END DO
      END DO
    END IF
  END FUNCTION

  SUBROUTINE MatrixGPUPlan(matrix,n,pr,pc,nnz_a,nnz_x,iterations,peak,seconds,safe)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    INTEGER, INTENT(IN) :: n,pr,pc,iterations
    INTEGER(INT64), INTENT(IN) :: nnz_a,nnz_x
    INTEGER(INT64), INTENT(OUT) :: peak
    REAL(REAL64), INTENT(OUT) :: seconds
    LOGICAL, INTENT(OUT) :: safe
    IF(matrix%is_complex) THEN
      CALL EstimateComplexProductPlan(n,pr,pc,nnz_a,nnz_x,iterations,peak,seconds,safe)
    ELSE
      CALL EstimateGPUPlan(n,pr,pc,nnz_a,nnz_x,iterations,peak,seconds,safe)
    END IF
  END SUBROUTINE

  SUBROUTINE ComplexCoordinates(matrix,root_mode,rows,columns,values,count)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    LOGICAL, INTENT(IN) :: root_mode
    INTEGER(C_INT), ALLOCATABLE, INTENT(OUT) :: rows(:),columns(:)
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, INTENT(OUT) :: values(:)
    INTEGER(C_INT), INTENT(OUT) :: count
    TYPE(TripletList_c) :: entries
    TYPE(Matrix_lsc) :: local
    INTEGER :: i
    CALL ConstructTripletList(entries)
    IF(root_mode) THEN
      CALL GatherMatrixToProcess(matrix,local,0,include_padding_in=.TRUE.)
      IF(matrix%process_grid%global_rank==0) CALL MatrixToTripletList(local,entries)
      CALL DestructMatrix(local)
    ELSE
      CALL GetMatrixTripletList(matrix,entries)
    END IF
    count=entries%CurrentSize
    ALLOCATE(rows(MAX(1,count)),columns(MAX(1,count)),values(MAX(1,count)))
    DO i=1,count
      rows(i)=entries%DATA(i)%index_row;columns(i)=entries%DATA(i)%index_column
      values(i)=entries%DATA(i)%point_value
    END DO
    CALL DestructTripletList(entries)
  END SUBROUTINE

  LOGICAL FUNCTION EvictCachedGPUWorkspaces(matrix) RESULT(evicted)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    INTEGER(C_INT) :: cleared,status
    INTEGER :: local_hit,global_hit,ierr
    status=trim_workspaces(cleared)
    IF(status/=0) CALL FailGPU('idle workspace reclamation failed')
    local_hit=INT(cleared)
    CALL MPI_Allreduce(local_hit,global_hit,1,MPI_INTEGER,MPI_MAX,matrix%process_grid%global_comm,ierr)
    evicted=global_hit/=0
  END FUNCTION

  LOGICAL FUNCTION EvictComplexWorkspace(matrix,limit) RESULT(evicted)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    INTEGER(C_INT64_T), INTENT(IN) :: limit
    INTEGER(C_INT) :: cleared,status
    INTEGER :: local_hit,global_hit,mode,ierr
    evicted=.FALSE.
    IF(limit<0) RETURN
    local_hit=0
    DO mode=1,2
      IF(.NOT.C_ASSOCIATED(matrix%process_grid%gpu_context(mode))) CYCLE
      status=clear_context_workspace(matrix%process_grid%gpu_context(mode),limit,cleared)
      IF(status/=0) CALL FailGPU('inactive workspace eviction failed')
      local_hit=MAX(local_hit,INT(cleared))
    END DO
    CALL MPI_Allreduce(local_hit,global_hit,1,MPI_INTEGER,MPI_MAX,matrix%process_grid%global_comm,ierr)
    evicted=global_hit/=0
  END FUNCTION

  SUBROUTINE TryComplexMultiply(a,b,c,alpha,beta,threshold,performed)
    TYPE(Matrix_ps), INTENT(IN) :: a,b
    TYPE(Matrix_ps), INTENT(INOUT) :: c
    REAL(NTREAL), INTENT(IN) :: alpha,beta,threshold
    LOGICAL, INTENT(OUT) :: performed
    TYPE(TripletList_c) :: entries
    INTEGER(C_INT), ALLOCATABLE :: ar(:),ac(:),br(:),bc(:),cr(:),cc(:),outrows(:),outcols(:)
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: av(:),bv(:),cv(:),outvalues(:)
    INTEGER(C_INT) :: na,nb,nc,count,rows,columns,n,row_offset,column_offset,status
    INTEGER :: valid,all_valid,ierr,i
    INTEGER(INT64) :: extra,global_extra,limit,trim_limit
    TYPE(C_PTR) :: output
    LOGICAL :: allowed,root_mode,participant,keep_c
    CHARACTER(16) :: report
    performed=.FALSE.;valid=1;keep_c=ABS(beta)>=TINY(beta)
    IF(nccl_active .OR. .NOT.ASSOCIATED(a%process_grid,b%process_grid)) valid=0
    IF(a%logical_matrix_dimension/=b%logical_matrix_dimension .OR. &
       a%actual_matrix_dimension/=b%actual_matrix_dimension) valid=0
    IF(.NOT.IEEE_IS_FINITE(alpha) .OR. .NOT.IEEE_IS_FINITE(beta) .OR. &
       .NOT.IEEE_IS_FINITE(threshold) .OR. threshold<0) valid=0
    IF(keep_c) THEN
      IF(.NOT.ASSOCIATED(a%process_grid,c%process_grid)) valid=0
      IF(a%logical_matrix_dimension/=c%logical_matrix_dimension .OR. &
         a%actual_matrix_dimension/=c%actual_matrix_dimension) valid=0
    END IF
    CALL MPI_Allreduce(valid,all_valid,1,MPI_INTEGER,MPI_MIN,a%process_grid%global_comm,ierr)
    IF(all_valid==0) RETURN
    ! Select using a complex input so the memory estimate includes complex payloads.
    extra=LocalGPUCount(b)
    IF(.NOT.a%is_complex) extra=LocalGPUCount(a)
    IF(keep_c) extra=extra+LocalGPUCount(c)
    CALL MPI_Allreduce(extra,global_extra,1,MPINTLONG,MPI_SUM,a%process_grid%global_comm,ierr)
    IF(a%is_complex) THEN
      CALL SelectGPUBackend(a,a,allowed,root_mode,1,global_extra,preserve_padding=.TRUE.,memory_limit=limit)
    ELSE
      CALL SelectGPUBackend(b,b,allowed,root_mode,1,global_extra,preserve_padding=.TRUE.,memory_limit=limit)
    END IF
    trim_limit=limit
    IF(.NOT.allowed) trim_limit=0
    IF(EvictComplexWorkspace(a,trim_limit)) THEN
      IF(a%is_complex) THEN
        CALL SelectGPUBackend(a,a,allowed,root_mode,1,global_extra,preserve_padding=.TRUE.)
      ELSE
        CALL SelectGPUBackend(b,b,allowed,root_mode,1,global_extra,preserve_padding=.TRUE.)
      END IF
    END IF
    IF(.NOT.allowed) RETURN
    CALL BeginGPUBackend(a,root_mode)
    IF(a%is_complex) THEN
      CALL CheckGPUMemory(a,a,root_mode,allowed,global_extra,preserve_padding=.TRUE.)
    ELSE
      CALL CheckGPUMemory(b,b,root_mode,allowed,global_extra,preserve_padding=.TRUE.)
    END IF
    participant=.NOT.root_mode .OR. a%process_grid%global_rank==0
    IF(.NOT.allowed) THEN
      IF(participant) THEN
        IF(clear_active_workspace()/=0) CALL FailGPU('active workspace eviction failed')
      END IF
      IF(a%is_complex) THEN
        CALL CheckGPUMemory(a,a,root_mode,allowed,global_extra,preserve_padding=.TRUE.)
      ELSE
        CALL CheckGPUMemory(b,b,root_mode,allowed,global_extra,preserve_padding=.TRUE.)
      END IF
    END IF
    IF(allowed) THEN
      CALL ComplexCoordinates(a,root_mode,ar,ac,av,na)
      CALL ComplexCoordinates(b,root_mode,br,bc,bv,nb)
      nc=0
      IF(keep_c) THEN
        CALL ComplexCoordinates(c,root_mode,cr,cc,cv,nc)
      ELSE
        ALLOCATE(cr(1),cc(1),cv(1))
      END IF
      CALL ConstructTripletList(entries)
      IF(participant) THEN
        n=a%logical_matrix_dimension;rows=a%local_rows;columns=a%local_columns
        row_offset=a%start_row-1;column_offset=a%start_column-1
        IF(root_mode) THEN
          rows=n;columns=n;row_offset=0;column_offset=0
        END IF
        status=complex_product(rows,columns,n,row_offset,column_offset,na,ar,ac,av,nb,br,bc,bv, &
          nc,cr,cc,cv,alpha,beta,threshold,output,count)
        IF(status/=0) CALL FailGPU('complex product failed')
        ALLOCATE(outrows(MAX(1,count)),outcols(MAX(1,count)),outvalues(MAX(1,count)))
        CALL complex_export(output,outrows,outcols,outvalues)
        CALL ConstructTripletList(entries,count)
        DO i=1,count
          entries%DATA(i)%index_row=outrows(i);entries%DATA(i)%index_column=outcols(i)
          entries%DATA(i)%point_value=outvalues(i)
        END DO
      END IF
      CALL ConstructEmptyMatrix(c,a%actual_matrix_dimension,process_grid_in=a%process_grid,is_complex_in=.TRUE.)
      CALL FillMatrixFromTripletList(c,entries,preduplicated_in=.TRUE.,prepartitioned_in=.NOT.root_mode)
      CALL DestructTripletList(entries)
      performed=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_ALGEBRA_REPORT',report)
      IF(TRIM(report)=='1' .AND. a%process_grid%global_rank==0) PRINT *, 'NTPoly GPU algebra complex product'
    END IF
    IF(participant) CALL EndGPUBackend()
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
    CALL ReleaseGPU(ga); CALL ReleaseGPU(gp); CALL ReleaseGPU(gpinv)
    CALL ReleaseGPU(temp); CALL ReleaseGPU(output)
    IF(.NOT.root_mode .OR. a%process_grid%global_rank==0) CALL EndGPUBackend()
  END SUBROUTINE
END MODULE
