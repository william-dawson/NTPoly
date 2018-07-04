!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Solve the matrix equation AX = B
MODULE LinearSolversModule
  USE CholeskyModule
  USE DataTypesModule
  USE MatrixDModule
  USE MatrixMemoryPoolPModule
  USE MatrixPSAlgebraModule
  USE MatrixPSModule
  USE FixedSolversModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE MatrixSModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CGSolver
  PUBLIC :: CholeskyDecomposition
  PUBLIC :: PivotedCholeskyDecomposition
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Solve the matrix equation AX = B using the conjugate gradient method.
  !! @param[in] AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] XMat the solved for matrix X.
  !! @param[in] BMat the right hand side.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE CGSolver(AMat, XMat, BMat, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: AMat
    TYPE(Matrix_ps), INTENT(INOUT) :: XMat
    TYPE(Matrix_ps), INTENT(IN)  :: BMat
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
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
       solver_parameters = IterativeSolverParameters_t()
    END IF

    !! Print out parameters
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Linear Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="CG")
       CALL PrintIterativeSolverParameters(solver_parameters)
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
    CALL StartTimer("Load Balance")
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
    CALL StopTimer("Load Balance")

    !! Initial Matrix Values
    CALL CopyMatrix(Identity, XMat)
    !! Compute residual
    CALL MatrixMultiply(ABalanced, Xmat, TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL CopyMatrix(BBalanced,RMat)
    CALL IncrementMatrix(TempMat, RMat, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL CopyMatrix(RMat,PMat)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
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
       CALL IncrementMatrix(QMat, RMat, alpha_in=-1.0*step_size)

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
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL PrintMatrixInformation(XMat)
    END IF

    !! Undo Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(XMat,XMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF
    CALL StopTimer("Load Balance")

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
  !> Compute The Cholesky Decomposition of a Symmetric Positive Definite matrix.
  !! This is a really naive implementation, that might be worth visiting.
  !! @param[in] AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] LMat the matrix computed.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE CholeskyDecomposition(AMat, LMat, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: AMat
    TYPE(Matrix_ps), INTENT(INOUT) :: LMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
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
       solver_parameters = FixedSolverParameters_t()
    END IF

    !! Print out parameters
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Linear Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Cholesky Decomposition")
       CALL PrintFixedSolverParameters(solver_parameters)
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
             inverse_factor = 1.0/insert_value
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
  !> Compute The Pivoted Cholesky Decomposition of a Symmetric Semi-Definite matrix.
  !! Pass it either a maximum rank and maximum error tolerance to check against
  !! @param[in] AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] LMat the matrix computed.
  !! @param[in] rank_in the target rank of the matrix.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE PivotedCholeskyDecomposition(AMat, LMat, rank_in, &
       & solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: AMat
    TYPE(Matrix_ps), INTENT(INOUT) :: LMat
    INTEGER, INTENT(IN) :: rank_in
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! For Pivoting
    INTEGER, DIMENSION(:), ALLOCATABLE :: pivot_vector
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: diag
    INTEGER :: pi_j
    !! Local Variables
    TYPE(Matrix_lsr) :: sparse_a
    TYPE(Matrix_lsr) :: acol
    TYPE(Matrix_ldr) :: dense_a
    !! For Storing The Results
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_column_l
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_l
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: values_l
    !! Broadcasted Column
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_index
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_values
    INTEGER :: recv_num_values
    !! Which pivots to treat locally
    INTEGER, DIMENSION(:), ALLOCATABLE :: local_pivots
    INTEGER :: num_local_pivots
    !! Root Lookups
    INTEGER, DIMENSION(:), ALLOCATABLE :: col_root_lookup
    !! Temporary Variables
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: a_buf
    INTEGER :: II, JJ, local_JJ
    INTEGER :: local_pi_i, local_pi_j
    REAL(NTREAL) :: Aval, insert_value, inverse_factor
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: dot_values
    INTEGER :: col_root
    INTEGER :: ierr

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF

    !! Print out parameters
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Linear Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", &
            & text_value_in="Pivoted Cholesky Decomposition")
       CALL WriteElement(key="Target_Rank", int_value_in=rank_in)
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    CALL ConstructEmptyMatrix(LMat, AMat)

    !! Construct the pivot vector
    ALLOCATE(pivot_vector(AMat%actual_matrix_dimension))
    DO JJ = 1, AMat%actual_matrix_dimension
       pivot_vector(JJ) = JJ
    END DO

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL MergeMatrixLocalBlocks(AMat, sparse_a)
    CALL ConstructMatrixDFromS(sparse_a, dense_a)
    ALLOCATE(local_pivots(sparse_a%columns))

    !! Extract the diagonal
    ALLOCATE(diag(sparse_a%columns))
    CALL ConstructDiag(AMat, Lmat%process_grid, dense_a, diag)

    !! Root Lookups
    ALLOCATE(col_root_lookup(AMat%logical_matrix_dimension))
    CALL ConstructRankLookup(AMat, LMat%process_grid, col_root_lookup)

    !! Allocate space for L
    ALLOCATE(values_per_column_l(sparse_a%columns))
    ALLOCATE(index_l(sparse_a%rows, sparse_a%columns))
    ALLOCATE(values_l(sparse_a%rows, sparse_a%columns))
    values_per_column_l = 0

    !! Buffer for fast indexing of A
    ALLOCATE(a_buf(sparse_a%columns))
    a_buf = 0

    !! Allocate space for a received column
    ALLOCATE(recv_index(sparse_a%rows))
    ALLOCATE(recv_values(sparse_a%rows))
    ALLOCATE(dot_values(sparse_a%columns))

    !! Pregather the full column of A.
    CALL GatherMatrixColumn(sparse_a, acol, LMat%process_grid)

    !! Main Loop
    DO JJ = 1, rank_in
       !! Pick a pivot vector
       local_JJ  = JJ - AMat%start_row + 1
       CALL GetPivot(AMat, LMat%process_grid, JJ, pivot_vector, diag, pi_j, &
            & insert_value, local_pivots, num_local_pivots)

       !! l[pi[j],j] = sqrt(d[pi[j]])
       IF (pi_j .GE. AMat%start_column .AND. pi_j .LT. AMat%end_column) THEN
          local_pi_j = pi_j - AMat%start_column + 1
          insert_value = SQRT(insert_value)
          inverse_factor = 1.0/insert_value
          !! Insert
          IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
             CALL AppendToVector(values_per_column_l(local_pi_j), &
                  & index_l(:,local_pi_j), values_l(:, local_pi_j), &
                  & local_JJ, insert_value)
          END IF
          !! Extract column J for sending later
          recv_num_values = values_per_column_l(local_pi_j)
          recv_index(:recv_num_values) = index_l(:recv_num_values,local_pi_j)
          recv_values(:recv_num_values) = values_l(:recv_num_values, local_pi_j)
       END IF
       col_root = col_root_lookup(pi_j)

       !! Broadcast column JJ, and Inverse Factor
       CALL BroadcastVector(recv_num_values, recv_index, recv_values, &
            & col_root, LMat%process_grid%row_comm)
       CALL MPI_Bcast(inverse_factor, 1, MPINTREAL, col_root, &
            & LMat%process_grid%row_comm, ierr)

       !! Extract the row of A to a dense matrix for easy lookup
       DO II = MAX(acol%outer_index(pi_j),1), acol%outer_index(pi_j+1)
          a_buf(acol%inner_index(II)) = acol%values(II)
       END DO

       !! Compute Dot Products
       CALL DotAllPivoted(recv_num_values, recv_index, recv_values, &
            & values_per_column_l, index_l, values_l, local_pivots, &
            & num_local_pivots, dot_values, LMat%process_grid%column_comm)

       !! Loop over other columns
       DO II = 1, num_local_pivots
          !! Insert Into L
          local_pi_i = local_pivots(II)
          Aval = a_buf(local_pi_i)
          insert_value = inverse_factor * (Aval - dot_values(II))
          IF (ABS(insert_value) .GT. solver_parameters%threshold) THEN
             IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
                CALL AppendToVector(values_per_column_l(local_pi_i), &
                     & index_l(:,local_pi_i), values_l(:, local_pi_i), &
                     & local_JJ, insert_value)
             END IF
          END IF
          !! Update Diagonal Array
          diag(local_pi_i) = diag(local_pi_i) - insert_value**2
       END DO

       !! Clear up the A buffer
       DO II = MAX(acol%outer_index(pi_j),1), acol%outer_index(pi_j+1)
          a_buf(acol%inner_index(II)) = 0
       END DO

    END DO

    !! Unpack
    CALL UnpackCholesky(values_per_column_l, index_l, values_l, LMat)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(LMat)
       CALL ExitSubLog
    END IF
    DEALLOCATE(local_pivots)
    DEALLOCATE(pivot_vector)
    DEALLOCATE(diag)
    CALL DestructMatrix(dense_a)
    DEALLOCATE(values_per_column_l)
    DEALLOCATE(index_l)
    DEALLOCATE(values_l)
    DEALLOCATE(recv_index)
    DEALLOCATE(recv_values)
    DEALLOCATE(dot_values)
    DEALLOCATE(a_buf)
    DEALLOCATE(col_root_lookup)
    CALL DestructMatrix(sparse_a)
  END SUBROUTINE PivotedCholeskyDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LinearSolversModule
