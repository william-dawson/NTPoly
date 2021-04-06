!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Solve the matrix equation AX = B
MODULE LinearSolversModule
  USE CholeskyModule, ONLY : ConstructRankLookup, AppendToVector, &
       & BroadcastVector, ConstructDiag, DotAllHelper, DotAllPivoted, &
       & GatherMatrixColumn, GetPivot, UnpackCholesky
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
       & FillMatrixIdentity, PrintMatrixInformation, MergeMatrixLocalBlocks
  USE SMatrixModule, ONLY : Matrix_lsr
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  USE NTMPIMODULE
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CGSolver
  PUBLIC :: CholeskyDecomposition
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
       CALL WriteElement(key="Method", VALUE="CG")
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
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
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
       CALL WriteElement(key="Total_Iterations", VALUE=outer_counter-1)
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
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE CGSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Hermitian Positive Definite matrix.
  !> This is a really naive implementation, that might be worth revisiting.
  SUBROUTINE CholeskyDecomposition(AMat, LMat, solver_parameters_in)
    !> The matrix A, must be hermitian, positive definite.
    TYPE(Matrix_ps), INTENT(IN)  :: AMat
    !> The lower diagonal matrix computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: LMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_lsr) :: sparse_a
    TYPE(Matrix_ldr) :: dense_a
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_column_l
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_l
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: values_l
    !! Temporary Variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_index
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_values
    INTEGER :: recv_num_values
    INTEGER :: II, JJ, local_JJ, local_II, local_row
    REAL(NTREAL) :: Aval, insert_value, inverse_factor
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: dot_values
    INTEGER, DIMENSION(:), ALLOCATABLE :: col_root_lookup
    INTEGER :: col_root
    INTEGER :: ierr

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
       CALL WriteElement(key="Method", value="Cholesky Decomposition")
       CALL PrintParameters(solver_parameters)
    END IF

    CALL ConstructEmptyMatrix(LMat, AMat)

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL MergeMatrixLocalBlocks(AMat, sparse_a)
    CALL ConstructMatrixDFromS(sparse_a, dense_a)

    !! Root Lookups
    ALLOCATE(col_root_lookup(AMat%logical_matrix_dimension))
    CALL ConstructRankLookup(AMat, LMat%process_grid, &
         & col_root_lookup=col_root_lookup)

    !! Allocate space for L
    ALLOCATE(values_per_column_l(sparse_a%columns))
    ALLOCATE(index_l(sparse_a%rows, sparse_a%columns))
    ALLOCATE(values_l(sparse_a%rows, sparse_a%columns))
    values_per_column_l = 0

    !! Allocate space for a received column
    ALLOCATE(recv_index(sparse_a%rows))
    ALLOCATE(recv_values(sparse_a%rows))
    ALLOCATE(dot_values(sparse_a%columns))

    !! Main Loop
    DO JJ = 1, AMat%actual_matrix_dimension
       !! Dot Column JJ with Column JJ, Insert Value into L[J,J]
       inverse_factor = 0
       IF (JJ .GE. AMat%start_column .AND. JJ .LT. AMat%end_column) THEN
          local_JJ = JJ - AMat%start_column + 1
          CALL DotAllHelper(values_per_column_l(local_JJ), &
               & index_l(:,local_JJ), values_l(:,local_JJ), &
               & values_per_column_l(local_JJ:local_JJ), &
               & index_l(:,local_JJ:local_JJ), values_l(:,local_JJ:local_JJ), &
               & dot_values(1:1), LMat%process_grid%column_comm)
          IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
             local_row = JJ - AMat%start_row + 1
             Aval = dense_a%data(local_row, local_JJ)
             insert_value = SQRT(Aval - dot_values(1))
             inverse_factor = 1.0_NTREAL/insert_value
             !! Insert
             CALL AppendToVector(values_per_column_l(local_JJ), &
                  & index_l(:,local_JJ), values_l(:, local_JJ), local_row, &
                  & insert_value)
          END IF
          !! Extract column J for sending later
          recv_num_values = values_per_column_l(local_JJ)
          recv_index(:recv_num_values) = index_l(:recv_num_values,local_JJ)
          recv_values(:recv_num_values) = values_l(:recv_num_values, local_JJ)
       END IF
       col_root = col_root_lookup(JJ)

       !! Broadcast column JJ, and Inverse Factor
       CALL BroadcastVector(recv_num_values, recv_index, recv_values, &
            & col_root, LMat%process_grid%row_comm)
       CALL MPI_Allreduce(MPI_IN_PLACE, inverse_factor, 1, MPINTREAL, MPI_SUM, &
            & LMat%process_grid%within_slice_comm, ierr)

       !! Loop over other columns
       CALL DotAllHelper(recv_num_values, recv_index, recv_values, &
            & values_per_column_l, index_l, values_l, dot_values, &
            & LMat%process_grid%column_comm)
       IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
          DO II = MAX(JJ + 1, AMat%start_column), AMat%end_column - 1
             local_II = II - AMat%start_column + 1
             local_row = JJ - AMat%start_row + 1
             Aval = dense_a%data(local_row, local_II)
             insert_value = inverse_factor * (Aval - dot_values(local_II))
             IF (ABS(insert_value) .GT. solver_parameters%threshold) THEN
                CALL AppendToVector(values_per_column_l(local_II), &
                     & index_l(:,local_II), values_l(:, local_II), &
                     & local_row, insert_value)
             END IF
          END DO
       END IF
    END DO

    !! Unpack
    CALL UnpackCholesky(values_per_column_l, index_l, values_l, LMat)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(LMat)
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(dense_a)
    DEALLOCATE(values_per_column_l)
    DEALLOCATE(index_l)
    DEALLOCATE(values_l)
    DEALLOCATE(recv_index)
    DEALLOCATE(recv_values)
    DEALLOCATE(dot_values)
    DEALLOCATE(col_root_lookup)
    CALL DestructMatrix(sparse_a)
  END SUBROUTINE CholeskyDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LinearSolversModule
