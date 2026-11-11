! An exact initial solution and scaled right-hand sides need finite CG results.
PROGRAM TestCGResidual
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE TripletModule
  USE TripletListModule
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE LinearSolversModule, ONLY: CGSolver
  USE SolverParametersModule
  IMPLICIT NONE
  TYPE(Matrix_ps) :: a,b,x,residual
  TYPE(SolverParameters_t) :: params
  TYPE(TripletList_r) :: entries
  TYPE(Triplet_r) :: entry
  REAL(NTREAL) :: factor,error
  INTEGER :: ierr,nr,rank,trial,j,k
  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr); CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  CALL ConstructProcessGrid(MPI_COMM_WORLD,1,nr,1)
  CALL ConstructEmptyMatrix(a,8); CALL FillMatrixIdentity(a); CALL ScaleMatrix(a,3.0_NTREAL)
  CALL ConstructSolverParameters(params,threshold_in=1e-8_NTREAL,converge_diff_in=1e-12_NTREAL, &
    max_iterations_in=40,monitor_convergence_in=.FALSE.)
  DO trial=1,4
    factor=1.0_NTREAL
    IF(trial==2) factor=1e-6_NTREAL
    IF(trial==3) factor=1e6_NTREAL
    IF(trial==4) THEN
      params%threshold=1e-12_NTREAL
      CALL ConstructTripletList(entries)
      IF(rank==0) THEN
        DO j=1,8
          entry%index_row=j; entry%index_column=j; entry%point_value=1.0_NTREAL+.1_NTREAL*j
          CALL AppendToTripletList(entries,entry)
          IF(j<8) THEN
            entry%index_column=j+1; entry%point_value=.0625_NTREAL
            CALL AppendToTripletList(entries,entry)
            entry%index_row=j+1; entry%index_column=j
            CALL AppendToTripletList(entries,entry)
          END IF
        END DO
      END IF
      CALL FillMatrixFromTripletList(a,entries); CALL DestructTripletList(entries)
    END IF
    CALL CopyMatrix(a,b); CALL ScaleMatrix(b,factor)
    IF(trial==4) CALL FillMatrixIdentity(b)
    CALL CGSolver(a,x,b,params)
    DO k=1,SIZE(x%local_data_r,2)
      DO j=1,SIZE(x%local_data_r,1)
        IF(ANY(.NOT.IEEE_IS_FINITE(x%local_data_r(j,k)%values))) CALL MPI_Abort(MPI_COMM_WORLD,2,ierr)
      END DO
    END DO
    CALL MatrixMultiply(a,x,residual,threshold_in=0.0_NTREAL)
    CALL IncrementMatrix(b,residual,alpha_in=-1.0_NTREAL)
    error=MatrixNorm(residual)/MatrixNorm(b)
    IF(.NOT.(error<1e-8_NTREAL)) THEN
      PRINT *, 'CG RESIDUAL FAILURE',trial,error
      CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF
    CALL DestructMatrix(b); CALL DestructMatrix(x); CALL DestructMatrix(residual)
  END DO
  CALL DestructMatrix(a); CALL DestructSolverParameters(params); CALL DestructProcessGrid()
  IF(rank==0) PRINT *, 'CG exact/scaled residual checks: PASS'
  CALL MPI_Finalize(ierr)
END PROGRAM
