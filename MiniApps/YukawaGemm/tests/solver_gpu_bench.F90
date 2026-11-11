! Complete public solver calls, including generic GPU import/export boundaries.
PROGRAM SolverGPUBench
  USE ISO_C_BINDING
  USE OMP_LIB
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE InverseSolversModule, ONLY: Invert
  USE SquareRootSolversModule, ONLY: InverseSquareRoot
  USE SolverParametersModule
  USE TripletListModule, ONLY: TripletList_r, DestructTripletList
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  TYPE(ProcessGrid_t), TARGET :: grids(0:2)
  TYPE(Matrix_ps) :: a(0:2),identity(0:2),reference,result,difference,temporary,residual
  TYPE(TripletList_r) :: reference_entries
  TYPE(SolverParameters_t) :: params
  CHARACTER(1024) :: filename,arg
  INTEGER :: ierr,rank,nr,pr,mode,solver,rep,status,provided,j,k
  LOGICAL :: gpu_only
  REAL(NTREAL) :: scale,error,t0,elapsed,residual_error
  INTERFACE
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_INT,C_CHAR
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
  END INTERFACE
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED,provided,ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr); CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,filename); CALL GET_COMMAND_ARGUMENT(2,arg); READ(arg,*) pr
  CALL GET_COMMAND_ARGUMENT(3,arg)
  gpu_only=TRIM(arg)=='gpu-only'
  ! Blocking is chosen at grid construction, so construct the CPU32 grid
  ! with32 threads rather than merely changing the thread count afterward.
  DO mode=0,2
    IF(gpu_only .AND. mode/=2) CYCLE
    CALL OMP_SET_NUM_THREADS(1)
    IF(mode==1) CALL OMP_SET_NUM_THREADS(32)
    CALL ConstructNewProcessGrid(grids(mode),MPI_COMM_WORLD,pr,nr/pr,1)
    CALL ConstructMatrixFromBinary(a(mode),TRIM(filename),process_grid_in=grids(mode))
    CALL ConstructEmptyMatrix(identity(mode),a(mode)); CALL FillMatrixIdentity(identity(mode))
    scale=MatrixNorm(a(mode))
    CALL ScaleMatrix(a(mode),0.2_NTREAL/scale); CALL IncrementMatrix(identity(mode),a(mode))
  END DO
  CALL ConstructSolverParameters(params,threshold_in=1e-8_NTREAL,converge_diff_in=1e-5_NTREAL, &
    max_iterations_in=50,monitor_convergence_in=.FALSE.)
  DO solver=1,2
    DO mode=0,2
      IF(gpu_only .AND. mode/=2) CYCLE
      status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
      CALL OMP_SET_NUM_THREADS(1)
      IF(mode==1) CALL OMP_SET_NUM_THREADS(32)
      IF(mode==2) status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,'1'//C_NULL_CHAR,1_C_INT)
      DO rep=1,3
        CALL MPI_Barrier(MPI_COMM_WORLD,ierr); t0=MPI_Wtime()
        IF(solver==1) CALL Invert(a(mode),result,params)
        IF(solver==2) CALL InverseSquareRoot(a(mode),result,params)
        elapsed=MPI_Wtime()-t0
        CALL MPI_Allreduce(MPI_IN_PLACE,elapsed,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        DO k=1,SIZE(result%local_data_r,2)
          DO j=1,SIZE(result%local_data_r,1)
            IF(ANY(.NOT.IEEE_IS_FINITE(result%local_data_r(j,k)%values))) CALL MPI_Abort(MPI_COMM_WORLD,4,ierr)
          END DO
        END DO
        ! Validate a mathematical residual outside the timed solver API.
        IF(solver==1) THEN
          CALL MatrixMultiply(a(mode),result,residual,threshold_in=0.0_NTREAL)
        ELSE
          CALL MatrixMultiply(a(mode),result,temporary,threshold_in=0.0_NTREAL)
          CALL MatrixMultiply(result,temporary,residual,threshold_in=0.0_NTREAL)
        END IF
        CALL IncrementMatrix(identity(mode),residual,alpha_in=-1.0_NTREAL)
        residual_error=MatrixNorm(residual)
        IF(rank==0) WRITE(*,'(A,3(1X,I0),1X,ES24.16)') 'solver_residual',solver,mode,rep,residual_error
        IF(.NOT.(residual_error<5e-5_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,5,ierr)
        CALL DestructMatrix(temporary); CALL DestructMatrix(residual)
        IF((mode==0 .OR. gpu_only) .AND. rep==1) CALL GetMatrixTripletList(result,reference_entries)
        CALL ConstructEmptyMatrix(reference,a(mode))
        CALL FillMatrixFromTripletList(reference,reference_entries)
        CALL CopyMatrix(result,difference)
        CALL IncrementMatrix(reference,difference,alpha_in=-1.0_NTREAL)
        error=MatrixNorm(difference)/MAX(1.0_NTREAL,MatrixNorm(reference))
        IF(rank==0) THEN
          WRITE(*,'(A,3(1X,I0),2(1X,ES24.16))') 'solver_bench',solver,mode,rep,elapsed,error
          FLUSH(6)
        END IF
        IF(.NOT.(error<1e-6_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,3,ierr)
        CALL DestructMatrix(result); CALL DestructMatrix(difference)
      END DO
    END DO
    CALL DestructMatrix(reference); CALL DestructTripletList(reference_entries)
  END DO
  DO mode=0,2
    IF(gpu_only .AND. mode/=2) CYCLE
    CALL DestructMatrix(a(mode)); CALL DestructMatrix(identity(mode)); CALL DestructProcessGrid(grids(mode))
  END DO
  CALL DestructSolverParameters(params); CALL MPI_Finalize(ierr)
END PROGRAM
