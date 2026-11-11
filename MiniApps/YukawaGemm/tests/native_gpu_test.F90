! Exercises the normal TRS2 API, owned matrices and host-return semantics.
PROGRAM NativeGPUTest
  USE ISO_C_BINDING
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule, ONLY: MatrixNorm, ScaleMatrix, IncrementMatrix
  USE DensityMatrixSolversModule, ONLY: TRS2
  USE SolverParametersModule
  USE PermutationModule
  USE TripletModule, ONLY: Triplet_r
  USE TripletListModule, ONLY: TripletList_r, ConstructTripletList, AppendToTripletList, DestructTripletList
  USE PSMatrixGPUBackendModule
  USE LoggingModule, ONLY: ActivateLogger, DeactivateLogger
  IMPLICIT NONE
  TYPE(Matrix_ps) :: h,identity,overlap,cpu,gpu,check_matrix
  TYPE(SolverParameters_t) :: params
  TYPE(Permutation_t) :: permutation
  TYPE(TripletList_r) :: entries
  TYPE(Triplet_r) :: entry
  CHARACTER(1024) :: filename,arg
  INTEGER :: ierr,rank,nr,pr,pc,example,j,repetition,status,work_comm
  REAL(NTREAL) :: ec,eg,muc,mug,error,scale
  INTERFACE
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_INT,C_CHAR
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
    INTEGER(C_INT) FUNCTION live_count() BIND(C,NAME='gpu_ops_live_count_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT) FUNCTION context_count() BIND(C,NAME='solver_nccl_context_count_')
      IMPORT C_INT
    END FUNCTION
  END INTERFACE
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr); CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,filename); CALL GET_COMMAND_ARGUMENT(2,arg); READ(arg,*) pr
  pc=nr/pr
  CALL GET_COMMAND_ARGUMENT(3,arg)
  IF(TRIM(arg)=='split') THEN
    CALL MPI_Comm_split(MPI_COMM_WORLD,MOD(rank,2),rank,work_comm,ierr)
    CALL MPI_Comm_rank(work_comm,rank,ierr); CALL MPI_Comm_size(work_comm,nr,ierr)
    pc=nr/pr
  ELSE
    CALL MPI_Comm_dup(MPI_COMM_WORLD,work_comm,ierr)
  END IF
  CALL ConstructProcessGrid(work_comm,pr,pc,1)
  CALL ConstructMatrixFromMatrixMarket(h,TRIM(filename))
  CALL ConstructEmptyMatrix(identity,h); CALL FillMatrixIdentity(identity)
  CALL ConstructDefaultPermutation(permutation,h%logical_matrix_dimension)
  CALL OwnershipTest()
  CALL BoundsTest()
  DO example=1,3
    CALL CopyMatrix(identity,overlap)
    IF(example==2) THEN
      CALL ConstructTripletList(entries)
      IF(rank==0) THEN
        DO j=1,h%actual_matrix_dimension
          entry%index_row=j; entry%index_column=j; entry%point_value=1.0_NTREAL+0.01_NTREAL*j
          CALL AppendToTripletList(entries,entry)
          IF(j<h%actual_matrix_dimension) THEN
            entry%index_row=j; entry%index_column=j+1; entry%point_value=0.005_NTREAL
            CALL AppendToTripletList(entries,entry)
            entry%index_row=j+1; entry%index_column=j
            CALL AppendToTripletList(entries,entry)
          END IF
        END DO
      END IF
      CALL FillMatrixFromTripletList(overlap,entries)
      CALL DestructTripletList(entries)
      ! Nontrivial deterministic permutation also exercises CPU setup/output fallback.
      DO j=1,h%actual_matrix_dimension
        permutation%index_lookup(j)=h%actual_matrix_dimension+1-j
        permutation%reverse_index_lookup(j)=h%actual_matrix_dimension+1-j
      END DO
    END IF
    IF(example==3) CALL ScaleMatrix(h,1.01_NTREAL)
    IF(example==2) THEN
      CALL ConstructSolverParameters(params,converge_diff_in=1e-10_NTREAL,threshold_in=1e-8_NTREAL, &
        max_iterations_in=60,be_verbose_in=.TRUE.,BalancePermutation_in=permutation)
    ELSE
      CALL ConstructSolverParameters(params,converge_diff_in=1e-10_NTREAL,threshold_in=1e-8_NTREAL, &
        max_iterations_in=60,be_verbose_in=.TRUE.)
    END IF
    DO repetition=1,2
      status=setenv('NTPOLY_GPU_TRS2'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
      IF(rank==0) CALL ActivateLogger()
      CALL TRS2(h,overlap,2.0_NTREAL,cpu,ec,muc,params)
      status=setenv('NTPOLY_GPU_TRS2'//C_NULL_CHAR,'1'//C_NULL_CHAR,1_C_INT)
      CALL TRS2(h,overlap,2.0_NTREAL,gpu,eg,mug,params)
      IF(rank==0) CALL DeactivateLogger()
      CALL CopyMatrix(gpu,check_matrix)
      CALL IncrementMatrix(cpu,check_matrix,alpha_in=-1.0_NTREAL)
      scale=MAX(1.0_NTREAL,MatrixNorm(cpu)); error=MatrixNorm(check_matrix)/scale
      IF(.NOT.(error<1e-11_NTREAL .AND. ABS(ec-eg)<1e-10_NTREAL .AND. ABS(muc-mug)<1e-9_NTREAL)) THEN
        PRINT *, 'NATIVE GPU PARITY FAILURE',rank,example,repetition,error,ec,eg,muc,mug
        CALL MPI_Abort(MPI_COMM_WORLD,7,ierr)
      END IF
      IF(live_count()/=0) CALL MPI_Abort(MPI_COMM_WORLD,8,ierr)
      IF(rank==0) PRINT *, 'NATIVE GPU PARITY PASS',example,repetition,error,ABS(ec-eg),ABS(muc-mug)
      CALL DestructMatrix(check_matrix); CALL DestructMatrix(cpu); CALL DestructMatrix(gpu)
    END DO
    CALL DestructSolverParameters(params); CALL DestructMatrix(overlap)
  END DO
  CALL ContextLifetimeTest()
  CALL DestructPermutation(permutation); CALL DestructMatrix(h); CALL DestructMatrix(identity)
  CALL DestructProcessGrid
  IF(context_count()/=0) CALL MPI_Abort(MPI_COMM_WORLD,11,ierr)
  IF(rank==0) PRINT *, 'GPU CONTEXT TEARDOWN PASS'
  CALL MPI_Comm_free(work_comm,ierr); CALL MPI_Finalize(ierr)
CONTAINS
  SUBROUTINE BoundsTest()
    TYPE(Matrix_ps) :: positive
    TYPE(GPUMatrix_t) :: device_matrix
    REAL(NTREAL) :: lower,upper,expected_lower
    LOGICAL :: root_mode
    INTEGER :: mode
    CALL ConstructEmptyMatrix(positive,5)
    CALL FillMatrixIdentity(positive)
    expected_lower=1.0_NTREAL
    IF(positive%logical_matrix_dimension>5) expected_lower=0.0_NTREAL
    DO mode=1,2
      root_mode=mode==2
      CALL BeginGPUBackend(positive,root_mode)
      CALL UploadGPURegion(positive,device_matrix,root_mode)
      CALL BoundsGPURegion(device_matrix,positive,root_mode,lower,upper)
      IF(lower/=expected_lower .OR. upper/=1.0_NTREAL) CALL MPI_Abort(MPI_COMM_WORLD,16,ierr)
      CALL ReleaseGPU(device_matrix)
      IF(.NOT.root_mode .OR. rank==0) CALL EndGPUBackend()
    END DO
    CALL DestructMatrix(positive)
    IF(rank==0) PRINT *, 'GPU BOUNDS COLUMN/PADDING PASS'
  END SUBROUTINE

  SUBROUTINE ContextLifetimeTest()
    TYPE(ProcessGrid_t), TARGET :: copied,reshaped
    TYPE(Matrix_ps) :: other
    INTEGER :: baseline
    baseline=context_count()
    ! Copy is an independent grid: no borrowed CUDA/NCCL ownership.
    CALL CopyProcessGrid(global_grid,copied)
    CALL ConstructEmptyMatrix(other,5,process_grid_in=copied)
    CALL BeginGPUBackend(other); CALL EndGPUBackend()
    IF(context_count()/=baseline+1) CALL MPI_Abort(MPI_COMM_WORLD,12,ierr)
    CALL BeginGPUBackend(h); CALL EndGPUBackend()
    CALL BeginGPUBackend(other); CALL EndGPUBackend()
    IF(context_count()/=baseline+1) CALL MPI_Abort(MPI_COMM_WORLD,13,ierr)
    CALL DestructMatrix(other); CALL DestructProcessGrid(copied)
    IF(context_count()/=baseline) CALL MPI_Abort(MPI_COMM_WORLD,14,ierr)
    ! Reuse the same Fortran object/MPI handle storage with a different grid.
    CALL ConstructNewProcessGrid(reshaped,work_comm,pc,pr,1)
    CALL ConstructEmptyMatrix(other,5,process_grid_in=reshaped)
    CALL BeginGPUBackend(other); CALL EndGPUBackend()
    CALL DestructMatrix(other); CALL DestructProcessGrid(reshaped)
    IF(context_count()/=baseline) CALL MPI_Abort(MPI_COMM_WORLD,15,ierr)
    IF(rank==0) PRINT *, 'GPU CONTEXT GRID COPY/SWITCH/RESHAPE PASS'
  END SUBROUTINE
  SUBROUTINE OwnershipTest()
    INTEGER :: i
    CALL BeginGPUBackend(h)
    BLOCK
      TYPE(GPUMatrix_t) :: owners(20),copies(20),one
      DO i=1,20
        CALL UploadGPU(h,owners(i))
      END DO
      copies=owners
      one=owners(1)
      CALL ScaleGPU(owners(1),2.0_NTREAL)
      CALL CopyMatrix(h,check_matrix); CALL DownloadGPU(one,check_matrix)
      CALL IncrementMatrix(h,check_matrix,alpha_in=-1.0_NTREAL)
      IF(MatrixNorm(check_matrix)>1e-14_NTREAL) CALL MPI_Abort(MPI_COMM_WORLD,9,ierr)
      CALL DestructMatrix(check_matrix)
    END BLOCK
    IF(live_count()/=0) CALL MPI_Abort(MPI_COMM_WORLD,10,ierr)
    CALL EndGPUBackend()
    IF(rank==0) PRINT *, 'OWNED GPU ARRAY ASSIGNMENT/FINALIZATION PASS'
  END SUBROUTINE
END PROGRAM
