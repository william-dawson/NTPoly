! Complete EigenDecomposition API timing, with spectrum fingerprints afterward.
PROGRAM EigenGPUBench
  USE ISO_C_BINDING
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE EigenSolversModule, ONLY: EigenDecomposition
  USE SolverParametersModule
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  TYPE(Matrix_ps) :: h,w,q,previous,difference
  TYPE(SolverParameters_t) :: params
  CHARACTER(1024) :: filename,arg
  INTEGER :: ierr,rank,nr,pr,reps,rep,mode,j,k
  REAL(NTREAL) :: t0,seconds,input_trace,input_square,trace,square,error
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,filename);CALL GET_COMMAND_ARGUMENT(2,arg);READ(arg,*) pr
  CALL GET_COMMAND_ARGUMENT(3,arg);READ(arg,*) reps
  CALL GET_COMMAND_ARGUMENT(4,arg);READ(arg,*) mode
  CALL ConstructProcessGrid(MPI_COMM_WORLD,pr,nr/pr,1)
  CALL ConstructMatrixFromBinary(h,TRIM(filename))
  CALL MatrixTrace(h,input_trace);CALL DotMatrix(h,h,input_square)
  CALL ConstructSolverParameters(params,threshold_in=1e-12_NTREAL,be_verbose_in=.FALSE.)
  DO rep=1,reps
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr);t0=MPI_Wtime()
    IF(mode==0) THEN
      CALL EigenDecomposition(h,w,solver_parameters_in=params)
    ELSE
      CALL EigenDecomposition(h,w,q,solver_parameters_in=params)
    END IF
    seconds=MPI_Wtime()-t0
    CALL MPI_Allreduce(MPI_IN_PLACE,seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    DO k=1,SIZE(w%local_data_r,2)
      DO j=1,SIZE(w%local_data_r,1)
        IF(ANY(.NOT.IEEE_IS_FINITE(w%local_data_r(j,k)%values))) CALL MPI_Abort(MPI_COMM_WORLD,3,ierr)
      END DO
    END DO
    CALL MatrixTrace(w,trace);CALL DotMatrix(w,w,square)
    error=MAX(ABS(trace-input_trace)/MAX(1.0_NTREAL,ABS(input_trace)), &
      ABS(square-input_square)/MAX(1.0_NTREAL,input_square))
    IF(.NOT.(error<1e-8_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,4,ierr)
    IF(rep==1) CALL CopyMatrix(w,previous)
    CALL CopyMatrix(w,difference);CALL IncrementMatrix(previous,difference,alpha_in=-1.0_NTREAL)
    error=MAX(error,MatrixNorm(difference)/MAX(1.0_NTREAL,MatrixNorm(previous)))
    IF(.NOT.(error<1e-8_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,5,ierr)
    IF(rank==0) THEN
      WRITE(*,'(A,2(1X,I0),4(1X,ES24.16))') 'eigen_bench',rep,mode,seconds,error,trace,square
      FLUSH(6)
    END IF
    CALL DestructMatrix(w);CALL DestructMatrix(q);CALL DestructMatrix(difference)
  END DO
  CALL DestructMatrix(previous);CALL DestructMatrix(h);CALL DestructSolverParameters(params)
  CALL DestructProcessGrid();CALL MPI_Finalize(ierr)
END PROGRAM
