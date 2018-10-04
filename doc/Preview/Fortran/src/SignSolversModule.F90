!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Matrix Sign Function.
MODULE SignSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteListElement, WriteElement, WriteCitation
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, &
       & MatrixNorm, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, CopyMatrix, DestructMatrix, &
       & FillMatrixIdentity, PrintMatrixInformation, TransposeMatrix, &
       & ConjugateMatrix, ConstructEmptyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SignFunction
  PUBLIC :: PolarDecomposition
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function.
  SUBROUTINE SignFunction(Mat, SignMat, solver_parameters_in)
    !> The input matrix.
    TYPE(Matrix_ps), INTENT(IN) :: Mat
    !> The sign of Mat.
    TYPE(Matrix_ps), INTENT(INOUT) :: SignMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Sign Function Solver")
       CALL EnterSubLog
       CALL WriteCitation("nicholas2008functions")
       CALL PrintParameters(solver_parameters)
    END IF

    CALL CoreComputation(Mat, SignMat, solver_parameters, .FALSE.)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE SignFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the polar decomposition of a matrix Mat = U*H.
  SUBROUTINE PolarDecomposition(Mat, Umat, Hmat, solver_parameters_in)
    !> The input matrix.
    TYPE(Matrix_ps), INTENT(IN) :: Mat
    !> The unitary polar factor.
    TYPE(Matrix_ps), INTENT(INOUT) :: Umat
    !> The hermitian matrix factor.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: Hmat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    TYPE(Matrix_ps) :: UmatT

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Polar Decomposition Solver")
       CALL EnterSubLog
       CALL WriteCitation("nicholas2008functions")
       CALL PrintParameters(solver_parameters)
    END IF

    CALL CoreComputation(Mat, Umat, solver_parameters, .TRUE.)

    IF (PRESENT(Hmat)) THEN
       CALL TransposeMatrix(Umat, UmatT)
       IF (UmatT%is_complex) THEN
          CALL ConjugateMatrix(UmatT)
       END IF
       CALL MatrixMultiply(UmatT, Mat, Hmat, &
            & threshold_in=solver_parameters%threshold)
       CALL DestructMatrix(UmatT)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE PolarDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This is the implementation routine for both the sign function and
  !> polar decomposition.
  SUBROUTINE CoreComputation(Mat, OutMat, solver_parameters, needs_transpose)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: Mat
    !> Output of the routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
    !> Whether we need to perform transposes in this routine (for polar).
    LOGICAL, INTENT(IN) :: needs_transpose
    !! Local Matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: Temp1
    TYPE(Matrix_ps) :: Temp2
    TYPE(Matrix_ps) :: OutMatT
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Data
    REAL(NTREAL), PARAMETER :: alpha = 1.69770248526_NTREAL
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: alpha_k
    REAL(NTREAL) :: xk
    REAL(NTREAL) :: norm_value
    INTEGER :: outer_counter

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Identity, Mat)
    CALL ConstructEmptyMatrix(Temp1, Mat)
    CALL ConstructEmptyMatrix(Temp2, Mat)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       !! Permute Matrices
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Mat, OutMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    ELSE
       CALL CopyMatrix(Mat,OutMat)
    END IF

    !! Initialize
    CALL GershgorinBounds(Mat,e_min,e_max)
    xk = ABS(e_min/e_max)
    CALL ScaleMatrix(OutMat,1.0_NTREAL/ABS(e_max))

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    iterate: DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       !! Update Scaling Factors
       alpha_k = MIN(SQRT(3.0_NTREAL/(1.0_NTREAL+xk+xk**2)), alpha)
       xk = 0.5_NTREAL*alpha_k*xk*(3.0_NTREAL-(alpha_k**2)*xk**2)

       IF (needs_transpose) THEN
          CALL TransposeMatrix(OutMat, OutMatT)
          IF (OutMatT%is_complex) THEN
             CALL ConjugateMatrix(OutMatT)
          END IF
          CALL MatrixMultiply(OutMatT, OutMat, Temp1, &
               & alpha_in=-1.0_NTREAL*alpha_k**2, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       ELSE
          CALL MatrixMultiply(OutMat, OutMat, Temp1, &
               & alpha_in=-1.0_NTREAL*alpha_k**2, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       END IF
       CALL IncrementMatrix(Identity,Temp1,alpha_in=3.0_NTREAL)

       CALL MatrixMultiply(OutMat, Temp1, Temp2, alpha_in=0.5_NTREAL*alpha_k, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

       CALL IncrementMatrix(Temp2, OutMat, alpha_in=-1.0_NTREAL)
       norm_value = MatrixNorm(OutMat)
       CALL CopyMatrix(Temp2,OutMat)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO iterate
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL PrintMatrixInformation(OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat,OutMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(OutMatT)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE CoreComputation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SignSolversModule
