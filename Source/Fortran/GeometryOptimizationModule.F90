!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Geometry Optimization
MODULE GeometryOptimizationModule
  USE DataTypesModule, ONLY : NTREAL
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteElement, WriteListElement
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, &
       & IncrementMatrix, ScaleMatrix, DotMatrix, SimilarityTransform, &
       & MatrixTrace
  USE PSMatrixModule, ONLY : Matrix_ps, DestructMatrix, ConstructEmptyMatrix, &
       & PrintMatrixInformation, CopyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  USE SquareRootSolversModule, ONLY : SquareRoot, InverseSquareRoot
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: PurificationExtrapolate
  PUBLIC :: LowdinExtrapolate
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a new guess at the Density Matrix after updating the geometry.
  !> Based on the purification algorithm in \cite niklasson2010trace .
  SUBROUTINE PurificationExtrapolate(PreviousDensity, Overlap, trace, &
       & NewDensity, solver_parameters_in)
    !> Previous density to extrapolate from.
    TYPE(Matrix_ps), INTENT(IN) :: PreviousDensity
    !> The overlap matrix of the new geometry.
    TYPE(Matrix_ps), INTENT(IN) :: Overlap
    !> The trace of the density matrix (usually the number of electrons).
    REAL(NTREAL), INTENT(IN) :: trace
    !> The extrapolated density.
    TYPE(Matrix_ps), INTENT(INOUT) :: NewDensity
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingDensity
    TYPE(Matrix_ps) :: WorkingOverlap
    TYPE(Matrix_ps) :: TempMat
    !! Local Variables
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Extrapolator")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Purification")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("niklasson2010trace")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(NewDensity, PreviousDensity)
    CALL ConstructEmptyMatrix(WorkingDensity, PreviousDensity)
    CALL ConstructEmptyMatrix(WorkingOverlap, PreviousDensity)

    !! Compute the working hamiltonian.
    CALL CopyMatrix(PreviousDensity, WorkingDensity)
    CALL CopyMatrix(Overlap, WorkingOverlap)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingDensity, WorkingDensity, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(WorkingOverlap, WorkingOverlap, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF

    II = 1
    DO II = 1, params%max_iterations
       !! Xn+1 = Xn S1 Xn
       CALL MatrixMultiply(WorkingDensity, WorkingOverlap, TempMat, &
            & threshold_in=params%threshold, memory_pool_in=pool)
       CALL MatrixMultiply(TempMat, WorkingDensity, NewDensity, &
            & threshold_in=params%threshold, memory_pool_in=pool)

       !! Figure Out Sigma Value
       CALL DotMatrix(WorkingDensity, WorkingOverlap, trace_value)

       !! Xn+1 = 2 Xn - Xn S1 Xn
       IF (trace .GT. trace_value) THEN
          CALL ScaleMatrix(NewDensity, -1.0_NTREAL)
          CALL IncrementMatrix(WorkingDensity, NewDensity, 2.0_NTREAL)
       END IF

       !! Check Convergence
       CALL IncrementMatrix(NewDensity, WorkingDensity, -1.0_NTREAL)
       norm_value = MatrixNorm(WorkingDensity)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
          CALL EnterSubLog
          CALL WriteElement(key="Trace", VALUE=trace_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF

       !! Xn = Xn+1
       CALL CopyMatrix(NewDensity, WorkingDensity)
    END DO
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", VALUE=II)
       CALL PrintMatrixInformation(NewDensity)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(NewDensity, NewDensity, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(WorkingDensity)
    CALL DestructMatrix(WorkingOverlap)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)

  END SUBROUTINE PurificationExtrapolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a new guess at the Density Matrix after updating the geometry.
  !> Based on the lowdin algorithm in \cite exner2002comparison .
  SUBROUTINE LowdinExtrapolate(PreviousDensity, OldOverlap, NewOverlap, &
       & NewDensity, solver_parameters_in)
    !> THe previous density to extrapolate from.
    TYPE(Matrix_ps), INTENT(IN)  :: PreviousDensity
    !> The old overlap matrix from the previous geometry.
    TYPE(Matrix_ps), INTENT(IN)  :: OldOverlap
    !> The new overlap matrix from the current geometry.
    TYPE(Matrix_ps), INTENT(IN)  :: NewOverlap
    !> The extrapolated density.
    TYPE(Matrix_ps), INTENT(INOUT) :: NewDensity
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: SQRMat
    TYPE(Matrix_ps) :: ISQMat
    TYPE(Matrix_ps) :: TempMat
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Extrapolator")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Lowdin")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("exner2002comparison")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    CALL SquareRoot(OldOverlap, SQRMat, params)
    CALL InverseSquareRoot(NewOverlap, ISQMat, params)

    CALL SimilarityTransform(PreviousDensity, SQRMat, SQRMat, TempMat, &
         & pool_in=pool, threshold_in=params%threshold)
    CALL SimilarityTransform(TempMat, ISQMat, ISQMat, NewDensity, &
         & pool_in=pool, threshold_in=params%threshold)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(SQRMat)
    CALL DestructMatrix(ISQMat)
    CALL DestructMatrix(TempMat)
    CALL DestructSolverParameters(params)

  END SUBROUTINE LowdinExtrapolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE GeometryOptimizationModule
