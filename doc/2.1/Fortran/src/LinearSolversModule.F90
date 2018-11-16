!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Solve the matrix equation AX = B
MODULE LinearSolversModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DMatrixModule, ONLY : Matrix_ldr, DestructMatrix, &
       & ConstructMatrixDFromS
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteHeader, WriteListElement
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix, MatrixNorm, &
       & MatrixMultiply, MatrixTrace, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & TransposeMatrix, DestructMatrix, ConjugateMatrix, CopyMatrix, &
       & FillMatrixIdentity, MergeMatrixLocalBlocks, PrintMatrixInformation
  USE SMatrixModule, ONLY : Matrix_lsr
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CGSolver
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Solve the matrix equation AX = B using the conjugate gradient method.
  SUBROUTINE CGSolver(AMat, XMat, BMat, solver_parameters_in)
    !> The matrix A, must be hermitian, positive definite.
    TYPE(Matrix_ps), INTENT(IN)  :: AMat
    !> The solved for matrix X.
    TYPE(Matrix_ps), INTENT(INOUT) :: XMat
    !> The right hand side.
    TYPE(Matrix_ps), INTENT(IN)  :: BMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: ABalanced
    TYPE(Matrix_ps) :: BBalanced
    TYPE(Matrix_ps) :: RMat, PMat, QMat
    TYPE(Matrix_ps) :: RMatT, PMatT
    TYPE(Matrix_ps) :: TempMat
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool
    REAL(NTREAL) :: top, bottom, new_top, step_size

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    !! Print out parameters
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Linear Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="CG")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Setup all the matrices
    CALL ConstructEmptyMatrix(Identity, AMat)
    CALL FillMatrixIdentity(Identity)
    CALL ConstructEmptyMatrix(ABalanced, AMat)
    CALL ConstructEmptyMatrix(BBalanced, AMat)
    CALL ConstructEmptyMatrix(RMat, AMat)
    CALL ConstructEmptyMatrix(PMat, AMat)
    CALL ConstructEmptyMatrix(QMat, AMat)
    CALL ConstructEmptyMatrix(TempMat, AMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(AMat, ABalanced, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BMat, BBalanced, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    ELSE
       CALL CopyMatrix(AMat,ABalanced)
       CALL CopyMatrix(BMat,BBalanced)
    END IF

    !! Initial Matrix Values
    CALL CopyMatrix(Identity, XMat)
    !! Compute residual
    CALL MatrixMultiply(ABalanced, Xmat, TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL CopyMatrix(BBalanced,RMat)
    CALL IncrementMatrix(TempMat, RMat, -1.0_NTREAL)
    CALL CopyMatrix(RMat,PMat)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", value=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", value=norm_value)
          CALL ExitSubLog
       END IF
       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF

       !! Compute the Step Size
       CALL MatrixMultiply(ABalanced, PMat, QMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

       CALL TransposeMatrix(RMat,RMatT)
       IF (RMatT%is_complex) THEN
          CALL ConjugateMatrix(RMatT)
       END IF
       CALL MatrixMultiply(RMatT, RMat, TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL MatrixTrace(TempMat, top)
       CALL TransposeMatrix(PMat,PMatT)
       IF (PMatT%is_complex) THEN
          CALL ConjugateMatrix(PMatT)
       END IF
       CALL MatrixMultiply(PMatT, QMat, TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL MatrixTrace(TempMat, bottom)
       step_size = top/bottom

       !! Update
       CALL IncrementMatrix(PMat, XMat, alpha_in=step_size)
       norm_value = ABS(step_size*MatrixNorm(PMat))
       CALL IncrementMatrix(QMat, RMat, alpha_in=-1.0_NTREAL*step_size)

       !! Update PMat
       CALL TransposeMatrix(RMat,RMatT)
       IF (RMatT%is_complex) THEN
          CALL ConjugateMatrix(RMatT)
       END IF
       CALL MatrixMultiply(RMatT, RMat, TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL MatrixTrace(TempMat, new_top)
       step_size = new_top / top
       CALL ScaleMatrix(PMat, step_size)
       CALL IncrementMatrix(RMat, PMat)

    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter-1)
       CALL PrintMatrixInformation(XMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(XMat,XMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(RMat)
    CALL DestructMatrix(PMat)
    CALL DestructMatrix(QMat)
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(ABalanced)
    CALL DestructMatrix(BBalanced)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE CGSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LinearSolversModule
