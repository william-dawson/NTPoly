!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Matrix Sign Function.
MODULE SignSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE EigenSolversModule, ONLY : DenseMatrixFunction
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteListElement, WriteElement
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, &
       & MatrixNorm, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, CopyMatrix, DestructMatrix, &
       & FillMatrixIdentity, PrintMatrixInformation, TransposeMatrix, &
       & ConjugateMatrix, ConstructEmptyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SignFunction
  PUBLIC :: DenseSignFunction
  PUBLIC :: PolarDecomposition
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function.
  SUBROUTINE SignFunction(InMat, OutMat, solver_parameters_in)
    !> The input matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InMat
    !> The sign of Mat.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Sign Function Solver")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("nicholas2008functions")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    CALL CoreComputation(InMat, OutMat, params, .FALSE.)

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
  END SUBROUTINE SignFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function (dense version).
  SUBROUTINE DenseSignFunction(InMat, OutputMat, solver_parameters_in)
    !> The matrix to compute the sign of.
    TYPE(Matrix_ps), INTENT(IN)  :: InMat
    !> The sign of the input matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Sign Function Solver")
       CALL EnterSubLog
    END IF

    !! Apply
    CALL DenseMatrixFunction(InMat, OutputMat, SignLambda, params)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE DenseSignFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the polar decomposition of a matrix Mat = U*H.
  SUBROUTINE PolarDecomposition(InMat, Umat, Hmat, solver_parameters_in)
    !> The input matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InMat
    !> The unitary polar factor.
    TYPE(Matrix_ps), INTENT(INOUT) :: Umat
    !> The hermitian matrix factor.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: Hmat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    TYPE(Matrix_ps) :: UmatT

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Polar Decomposition Solver")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("nicholas2008functions")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    CALL CoreComputation(InMat, Umat, params, .TRUE.)

    IF (PRESENT(Hmat)) THEN
       CALL TransposeMatrix(Umat, UmatT)
       IF (UmatT%is_complex) THEN
          CALL ConjugateMatrix(UmatT)
       END IF
       CALL MatrixMultiply(UmatT, InMat, Hmat, &
            & threshold_in=params%threshold)
       CALL DestructMatrix(UmatT)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
  END SUBROUTINE PolarDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This is the implementation routine for both the sign function and
  !> polar decomposition.
  SUBROUTINE CoreComputation(InMat, OutMat, params, needs_transpose)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: InMat
    !> Output of the routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: params
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
    INTEGER :: II

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Identity, InMat)
    CALL ConstructEmptyMatrix(Temp1, InMat)
    CALL ConstructEmptyMatrix(Temp2, InMat)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       !! Permute Matrices
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(InMat, OutMat, &
            & params%BalancePermutation, memorypool_in=pool)
    ELSE
       CALL CopyMatrix(InMat, OutMat)
    END IF

    !! Initialize
    CALL GershgorinBounds(InMat, e_min, e_max)
    xk = ABS(e_min/e_max)
    CALL ScaleMatrix(OutMat, 1.0_NTREAL/ABS(e_max))

    !! Iterate.
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    iterate: DO II = 1, params%max_iterations
       IF (params%be_verbose .AND. II .GT. 1) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
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
               & threshold_in=params%threshold, memory_pool_in=pool)
       ELSE
          CALL MatrixMultiply(OutMat, OutMat, Temp1, &
               & alpha_in=-1.0_NTREAL*alpha_k**2, &
               & threshold_in=params%threshold, memory_pool_in=pool)
       END IF
       CALL IncrementMatrix(Identity, Temp1, alpha_in=3.0_NTREAL)

       CALL MatrixMultiply(OutMat, Temp1, Temp2, alpha_in=0.5_NTREAL*alpha_k, &
            & threshold_in=params%threshold, memory_pool_in=pool)

       CALL IncrementMatrix(Temp2, OutMat, alpha_in=-1.0_NTREAL)
       norm_value = MatrixNorm(OutMat)
       CALL CopyMatrix(Temp2, OutMat)

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF
    END DO iterate
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total Iterations", VALUE=II-1)
       CALL PrintMatrixInformation(OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat,OutMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(OutMatT)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE CoreComputation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prototypical sign function for mapping. 
  FUNCTION SignLambda(val) RESULT(outval)
    REAL(KIND=NTREAL), INTENT(IN) :: val
    REAL(KIND=NTREAL) :: outval

    IF (val < 0.0_NTREAL) THEN
       outval = -1.0_NTREAL
    ELSE
       outval = 1.0_NTREAL
    END IF
  END FUNCTION SignLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SignSolversModule
