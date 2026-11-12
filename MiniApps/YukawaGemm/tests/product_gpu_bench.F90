! Matched real/complex public multiplication costs. All validation is untimed.
PROGRAM ProductGPUBench
  USE ISO_C_BINDING
  USE OMP_LIB
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL,MPINTLONG
  USE ProcessGridModule
  USE PSMatrixModule
  USE PMatrixMemoryPoolModule
  USE PSMatrixGPUBackendModule, ONLY: GetGPUProductTimings
  USE PSMatrixAlgebraModule
  USE TripletListModule
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  TYPE(ProcessGrid_t), TARGET :: reference_grid,grids(0:2)
  TYPE(MatrixMemoryPool_p) :: workspace_pool(0:2)
  TYPE(Matrix_ps) :: input,reference(0:2),a(0:2),result,difference
  TYPE(TripletList_r) :: input_entries,result_entries
  TYPE(TripletList_c) :: complex_entries
  CHARACTER(1024) :: filename,arg,backend
  INTEGER :: ierr,rank,nr,pr,n,reps,threads,rep,kind,step,status,provided,j,k,selected_kind
  LOGICAL :: profile,stage_timing,use_pool,check_workspace
  INTEGER(C_INT64_T) :: allocation_before
  REAL(C_DOUBLE) :: stage_local(43)
  REAL(C_DOUBLE), ALLOCATABLE :: stage_all(:,:)
  REAL(NTREAL) :: t0,seconds,error,threshold
  INTEGER(KIND=8) :: nnz_input,nnz_output
  INTERFACE
    SUBROUTINE profile_begin() BIND(C,NAME='gpu_profile_begin_')
    END SUBROUTINE
    SUBROUTINE profile_end() BIND(C,NAME='gpu_profile_end_')
    END SUBROUTINE
    INTEGER(C_INT64_T) FUNCTION workspace_allocations() BIND(C,NAME='gpu_workspace_allocation_count_')
      IMPORT C_INT64_T
    END FUNCTION
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_INT,C_CHAR
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
  END INTERFACE
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED,provided,ierr)
  IF(provided<MPI_THREAD_SERIALIZED) CALL MPI_Abort(MPI_COMM_WORLD,9,ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,filename);CALL GET_COMMAND_ARGUMENT(2,arg);READ(arg,*) pr
  CALL GET_COMMAND_ARGUMENT(3,arg);READ(arg,*) n
  CALL GET_COMMAND_ARGUMENT(4,backend)
  CALL GET_COMMAND_ARGUMENT(5,arg);READ(arg,*) threads
  CALL GET_COMMAND_ARGUMENT(6,arg);READ(arg,*) reps
  CALL GET_COMMAND_ARGUMENT(7,arg);profile=TRIM(arg)=='profile'
  selected_kind=-1
  CALL GET_COMMAND_ARGUMENT(8,arg)
  IF(LEN_TRIM(arg)>0) READ(arg,*) selected_kind
  IF(selected_kind< -1 .OR. selected_kind>2) CALL MPI_Abort(MPI_COMM_WORLD,8,ierr)
  CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_ALGEBRA_TIMING',arg)
  stage_timing=TRIM(arg)=='1' .AND. TRIM(backend)=='gpu'
  IF(stage_timing) ALLOCATE(stage_all(43,nr))
  CALL GET_COMMAND_ARGUMENT(9,arg);use_pool=TRIM(arg)=='pool'
  CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_WORKSPACE_CHECK',arg);check_workspace=TRIM(arg)=='1'
  threshold=1e-8_NTREAL
  status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
  CALL OMP_SET_NUM_THREADS(1)
  CALL ConstructNewProcessGrid(reference_grid,MPI_COMM_WORLD,pr,nr/pr,1)
  CALL ConstructMatrixFromBinary(input,TRIM(filename),process_grid_in=reference_grid)
  IF(n>0 .AND. n/=input%actual_matrix_dimension) CALL ResizeMatrix(input,n)
  n=input%actual_matrix_dimension
  CALL GetMatrixTripletList(input,input_entries)
  nnz_input=GetMatrixSize(input)
  CALL MatrixMultiply(input,input,reference(0),threshold_in=threshold)
  CALL GetMatrixTripletList(reference(0),result_entries)
  DO kind=1,2
    CALL PhaseEntries(result_entries,complex_entries,kind==2)
    CALL ConstructEmptyMatrix(reference(kind),n,process_grid_in=reference_grid,is_complex_in=.TRUE.)
    CALL FillMatrixFromTripletList(reference(kind),complex_entries)
    CALL DestructTripletList(complex_entries)
  END DO
  CALL DestructTripletList(result_entries)
  DO kind=0,2
    CALL OMP_SET_NUM_THREADS(threads)
    CALL ConstructNewProcessGrid(grids(kind),MPI_COMM_WORLD,pr,nr/pr,1)
    CALL ConstructEmptyMatrix(a(kind),n,process_grid_in=grids(kind),is_complex_in=kind/=0)
    IF(kind==0) THEN
      CALL FillMatrixFromTripletList(a(kind),input_entries)
    ELSE
      CALL PhaseEntries(input_entries,complex_entries,kind==2)
      CALL FillMatrixFromTripletList(a(kind),complex_entries)
      CALL DestructTripletList(complex_entries)
    END IF
  END DO
  IF(TRIM(backend)=='gpu') status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,'1'//C_NULL_CHAR,1_C_INT)
  IF(rank==0) THEN
    WRITE(*,'(A,4(1X,I0),1X,A,1X,ES24.16)') 'product_input',n,nr,threads,nnz_input,TRIM(backend),threshold
    FLUSH(6)
  END IF
  DO rep=1,reps
    IF(profile .AND. rep==2) CALL profile_begin()
    DO step=0,2
      ! Rotate ordering so one datatype does not always benefit from running last.
      kind=MOD(step+rep-1,3)
      IF(selected_kind>=0 .AND. selected_kind/=kind) CYCLE
      IF(check_workspace .AND. kind==0) allocation_before=workspace_allocations()
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr);t0=MPI_Wtime()
      IF(use_pool) THEN
        CALL MatrixMultiply(a(kind),a(kind),result,threshold_in=threshold,memory_pool_in=workspace_pool(kind))
      ELSE
        CALL MatrixMultiply(a(kind),a(kind),result,threshold_in=threshold)
      END IF
      seconds=MPI_Wtime()-t0
      IF(check_workspace .AND. kind==0 .AND. rep>1) THEN
        IF(workspace_allocations()/=allocation_before) CALL MPI_Abort(MPI_COMM_WORLD,11,ierr)
      END IF
      IF(stage_timing .AND. kind==0) THEN
        stage_local(1)=seconds
        CALL GetGPUProductTimings(stage_local(2:43))
        CALL MPI_Gather(stage_local,43,MPI_DOUBLE_PRECISION,stage_all,43,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        IF(rank==0) THEN
          DO j=1,nr
            WRITE(*,'(A,3(1X,I0),43(1X,ES24.16))') 'product_stages',kind,rep,j-1,stage_all(:,j)
          END DO
        END IF
      END IF
      CALL MPI_Allreduce(MPI_IN_PLACE,seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      CALL CheckFinite(result)
      CALL OMP_SET_NUM_THREADS(1)
      ! Move the result onto the reference grid before comparing. Different
      ! OpenMP block decompositions must not be mixed in IncrementMatrix.
      IF(kind==0) THEN
        CALL GetMatrixTripletList(result,result_entries)
        CALL ConstructEmptyMatrix(difference,n,process_grid_in=reference_grid)
        CALL FillMatrixFromTripletList(difference,result_entries)
        CALL DestructTripletList(result_entries)
      ELSE
        CALL GetMatrixTripletList(result,complex_entries)
        CALL ConstructEmptyMatrix(difference,n,process_grid_in=reference_grid,is_complex_in=.TRUE.)
        CALL FillMatrixFromTripletList(difference,complex_entries)
        CALL DestructTripletList(complex_entries)
      END IF
      nnz_output=GetMatrixSize(result)
      CALL IncrementMatrix(reference(kind),difference,alpha_in=-1.0_NTREAL)
      error=MatrixNorm(difference)/MAX(1.0_NTREAL,MatrixNorm(reference(kind)))
      IF(.NOT.(error<1e-9_NTREAL)) THEN
        IF(rank==0) PRINT *, 'PRODUCT BENCH VALIDATION FAILURE',kind,rep,error
        CALL MPI_Abort(MPI_COMM_WORLD,3,ierr)
      END IF
      IF(rank==0) THEN
        WRITE(*,'(A,2(1X,I0),2(1X,ES24.16),1X,I0)') 'product_bench',kind,rep,seconds,error,nnz_output
        FLUSH(6)
      END IF
      CALL DestructMatrix(result);CALL DestructMatrix(difference)
      CALL OMP_SET_NUM_THREADS(threads)
    END DO
  END DO
  IF(profile) CALL profile_end()
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr);t0=MPI_Wtime()
  DO kind=0,2
    CALL DestructMatrixMemoryPool(workspace_pool(kind))
    CALL DestructMatrix(a(kind));CALL DestructMatrix(reference(kind));CALL DestructProcessGrid(grids(kind))
  END DO
  CALL DestructTripletList(input_entries);CALL DestructMatrix(input);CALL DestructProcessGrid(reference_grid)
  seconds=MPI_Wtime()-t0
  CALL MPI_Allreduce(MPI_IN_PLACE,seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  IF(rank==0) THEN
    WRITE(*,'(A,1X,ES24.16)') 'product_teardown',seconds
    PRINT *, 'PRODUCT BENCH COMPLETE'
  END IF
  CALL MPI_Finalize(ierr)
CONTAINS
  SUBROUTINE PhaseEntries(real_entries,entries,phase)
    TYPE(TripletList_r), INTENT(IN) :: real_entries
    TYPE(TripletList_c), INTENT(INOUT) :: entries
    LOGICAL, INTENT(IN) :: phase
    INTEGER :: i,r,c
    REAL(NTREAL) :: theta
    CALL CopyTripletList(real_entries,entries)
    IF(.NOT.phase) RETURN
    DO i=1,entries%CurrentSize
      r=entries%DATA(i)%index_row;c=entries%DATA(i)%index_column
      theta=2.0_NTREAL*ACOS(-1.0_NTREAL)*REAL(MOD(97*r,257)-MOD(97*c,257),NTREAL)/257
      entries%DATA(i)%point_value=entries%DATA(i)%point_value*CMPLX(COS(theta),SIN(theta),NTREAL)
    END DO
  END SUBROUTINE
  SUBROUTINE CheckFinite(matrix)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    IF(matrix%is_complex) THEN
      DO k=1,SIZE(matrix%local_data_c,2)
        DO j=1,SIZE(matrix%local_data_c,1)
          IF(ANY(.NOT.IEEE_IS_FINITE(REAL(matrix%local_data_c(j,k)%values,NTREAL))) .OR. &
             ANY(.NOT.IEEE_IS_FINITE(AIMAG(matrix%local_data_c(j,k)%values)))) CALL MPI_Abort(MPI_COMM_WORLD,4,ierr)
        END DO
      END DO
    ELSE
      DO k=1,SIZE(matrix%local_data_r,2)
        DO j=1,SIZE(matrix%local_data_r,1)
          IF(ANY(.NOT.IEEE_IS_FINITE(matrix%local_data_r(j,k)%values))) CALL MPI_Abort(MPI_COMM_WORLD,4,ierr)
        END DO
      END DO
    END IF
  END SUBROUTINE
END PROGRAM
