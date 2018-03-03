!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Solve the matrix equation AX = B
MODULE LinearSolversModule
  USE DataTypesModule
  USE DenseMatrixModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixAlgebraModule
  USE DistributedSparseMatrixModule
  USE FixedSolversModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE MatrixGatherModule
  USE ProcessGridModule
  USE SparseMatrixAlgebraModule
  USE SparseMatrixModule
  USE SparseVectorModule
  USE TimerModule
  USE TripletListModule
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
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
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: AMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: XMat
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: BMat
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: ABalanced
    TYPE(DistributedSparseMatrix_t) :: BBalanced
    TYPE(DistributedSparseMatrix_t) :: RMat, PMat, QMat
    TYPE(DistributedSparseMatrix_t) :: RMatT, PMatT
    TYPE(DistributedSparseMatrix_t) :: TempMat
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(DistributedMatrixMemoryPool_t) :: pool
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
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & AMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)
    CALL ConstructEmptyDistributedSparseMatrix(ABalanced, &
         & AMat%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(BBalanced, &
         & AMat%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(RMat, &
         & AMat%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(PMat, &
         & AMat%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(QMat, &
         & AMat%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat, &
         & AMat%actual_matrix_dimension)

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
       CALL CopyDistributedSparseMatrix(AMat,ABalanced)
       CALL CopyDistributedSparseMatrix(BMat,BBalanced)
    END IF
    CALL StopTimer("Load Balance")

    !! Initial Matrix Values
    CALL CopyDistributedSparseMatrix(Identity, XMat)
    !! Compute residual
    CALL DistributedGemm(ABalanced, Xmat, TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL CopyDistributedSparseMatrix(BBalanced,RMat)
    CALL IncrementDistributedSparseMatrix(TempMat, RMat, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL CopyDistributedSparseMatrix(RMat,PMat)

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
       CALL DistributedGemm(ABalanced, PMat, QMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       !  top = DotDistributedSparseMatrix(RMat,RMat)
       !  bottom = DotDistributedSparseMatrix(PMat,QMat)
       !  step_size = top/bottom
       CALL TransposeDistributedSparseMatrix(RMat,RMatT)
       CALL DistributedGemm(RMatT, RMat, TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       top = Trace(TempMat)
       CALL TransposeDistributedSparseMatrix(PMat,PMatT)
       CALL DistributedGemm(PMatT, QMat, TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       bottom = Trace(TempMat)
       step_size = top/bottom

       !! Update
       CALL IncrementDistributedSparseMatrix(PMat, XMat, alpha_in=step_size)
       norm_value = ABS(step_size*DistributedSparseNorm(PMat))
       CALL IncrementDistributedSparseMatrix(QMat, RMat, alpha_in=-1.0*step_size)

       !! Update PMat
       !!new_top = DotDistributedSparseMatrix(RMat,RMat)
       CALL TransposeDistributedSparseMatrix(RMat,RMatT)
       CALL DistributedGemm(RMatT, RMat, TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       new_top = Trace(TempMat)
       step_size = new_top / top
       CALL ScaleDistributedSparseMatrix(PMat, step_size)
       CALL IncrementDistributedSparseMatrix(RMat, PMat)

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
    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL ExitSubLog
    END IF
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(RMat)
    CALL DestructDistributedSparseMatrix(PMat)
    CALL DestructDistributedSparseMatrix(QMat)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(ABalanced)
    CALL DestructDistributedSparseMatrix(BBalanced)
  END SUBROUTINE CGSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Symmetric Positive Definite matrix.
  !! This is a really naive implementation, that might be worth visiting.
  !! @param[in] AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] LMat the matrix computed.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE CholeskyDecomposition(AMat, LMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: AMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: LMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(SparseMatrix_t) :: sparse_a
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: dense_a
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

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL MergeLocalBlocks(AMat, sparse_a)
    ALLOCATE(dense_a(sparse_a%rows, sparse_a%columns))
    dense_a = 0
    CALL ConstructDenseFromSparse(sparse_a, dense_a)

    !! Root Lookups
    ALLOCATE(col_root_lookup(AMat%logical_matrix_dimension))
    CALL ConstructRankLookup(AMat, col_root_lookup=col_root_lookup)

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
               & dot_values(1:1), column_comm)
          IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
             local_row = JJ - AMat%start_row + 1
             Aval = dense_a(local_row, local_JJ)
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
            & col_root, row_comm)
       CALL MPI_Allreduce(MPI_IN_PLACE, inverse_factor, 1, MPINTREAL, MPI_SUM, &
            & within_slice_comm, grid_error)

       !! Loop over other columns
       CALL DotAllHelper(recv_num_values, recv_index, recv_values, &
            & values_per_column_l, index_l, values_l, dot_values, column_comm)
       IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
          DO II = MAX(JJ + 1, AMat%start_column), AMat%end_column - 1
             local_II = II - AMat%start_column + 1
             local_row = JJ - AMat%start_row + 1
             Aval = dense_a(local_row, local_II)
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
    CALL ConstructEmptyDistributedSparseMatrix(LMat, &
         & AMat%actual_matrix_dimension)
    CALL UnpackCholesky(values_per_column_l, index_l, values_l, LMat)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(LMat)
       CALL ExitSubLog
    END IF
    DEALLOCATE(dense_a)
    DEALLOCATE(values_per_column_l)
    DEALLOCATE(index_l)
    DEALLOCATE(values_l)
    DEALLOCATE(recv_index)
    DEALLOCATE(recv_values)
    DEALLOCATE(dot_values)
    DEALLOCATE(col_root_lookup)
    CALL DestructSparseMatrix(sparse_a)
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
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: AMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: LMat
    INTEGER, INTENT(IN) :: rank_in
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! For Pivoting
    INTEGER, DIMENSION(:), ALLOCATABLE :: pivot_vector
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: diag
    INTEGER :: pi_j
    !! Local Variables
    TYPE(SparseMatrix_t) :: sparse_a
    TYPE(SparseMatrix_t) :: acol
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: dense_a
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

    !! Construct the pivot vector
    ALLOCATE(pivot_vector(AMat%actual_matrix_dimension))
    DO JJ = 1, AMat%actual_matrix_dimension
       pivot_vector(JJ) = JJ
    END DO

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL MergeLocalBlocks(AMat, sparse_a)
    ALLOCATE(dense_a(sparse_a%rows, sparse_a%columns))
    dense_a = 0
    CALL ConstructDenseFromSparse(sparse_a, dense_a)
    ALLOCATE(local_pivots(sparse_a%columns))

    !! Extract the diagonal
    ALLOCATE(diag(sparse_a%columns))
    CALL ConstructDiag(AMat, dense_a, diag)

    !! Root Lookups
    ALLOCATE(col_root_lookup(AMat%logical_matrix_dimension))
    CALL ConstructRankLookup(AMat, col_root_lookup)

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
    CALL GatherMatrixColumn(sparse_a, acol)

    !! Main Loop
    DO JJ = 1, rank_in
       !! Pick a pivot vector
       local_JJ  = JJ - AMat%start_row + 1
       CALL GetPivot(AMat, JJ, pivot_vector, diag, pi_j, local_pivots, &
            & num_local_pivots)

       !! l[pi[j],j] = sqrt(d[pi[j]])
       IF (pi_j .GE. AMat%start_column .AND. pi_j .LT. AMat%end_column) THEN
          local_pi_j = pi_j - AMat%start_column + 1
          insert_value = SQRT(diag(local_pi_j))
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
            & col_root, row_comm)
       CALL MPI_Bcast(inverse_factor, 1, MPINTREAL, col_root, row_comm, &
            & grid_error)

       !! Extract the row of A to a dense matrix for easy lookup
       DO II = MAX(acol%outer_index(pi_j),1), acol%outer_index(pi_j+1)
          a_buf(acol%inner_index(II)) = acol%values(II)
       END DO

       !! Compute Dot Products
       CALL DotAllPivoted(recv_num_values, recv_index, recv_values, &
            & values_per_column_l, index_l, values_l, local_pivots, &
            & num_local_pivots, dot_values, column_comm)

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
    CALL ConstructEmptyDistributedSparseMatrix(LMat, &
         & AMat%actual_matrix_dimension)
    CALL UnpackCholesky(values_per_column_l, index_l, values_l, LMat)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(LMat)
       CALL ExitSubLog
    END IF
    DEALLOCATE(local_pivots)
    DEALLOCATE(pivot_vector)
    DEALLOCATE(diag)
    DEALLOCATE(dense_a)
    DEALLOCATE(values_per_column_l)
    DEALLOCATE(index_l)
    DEALLOCATE(values_l)
    DEALLOCATE(recv_index)
    DEALLOCATE(recv_values)
    DEALLOCATE(dot_values)
    DEALLOCATE(a_buf)
    DEALLOCATE(col_root_lookup)
    CALL DestructSparseMatrix(sparse_a)
  END SUBROUTINE PivotedCholeskyDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE UnpackCholesky(values_per_column, index, values, LMat)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN) :: values_per_column
    INTEGER, DIMENSION(:,:), INTENT(IN) :: index
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: LMat
    !! Local Variables
    INTEGER :: local_columns
    TYPE(TripletList_t) :: local_triplets
    TYPE(Triplet_t) :: temp
    INTEGER :: II, JJ

    local_columns = LMat%local_columns

    CALL ConstructTripletList(local_triplets)
    IF (my_slice .EQ. 0) THEN
       DO JJ = 1, local_columns
          !! note transpose
          temp%index_row = JJ + LMat%start_column - 1
          DO II = 1, values_per_column(JJ)
             !! note transpose
             temp%index_column = INDEX(II,JJ) + LMat%start_row - 1
             temp%point_value = values(II,JJ)
             CALL AppendToTripletList(local_triplets, temp)
          END DO
       END DO
    END IF
    CALL FillFromTripletList(LMat, local_triplets)

    !! Cleanup
    CALL DestructTripletList(local_triplets)
  END SUBROUTINE UnpackCholesky
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine which computes sparse dot products across processors.
  !! Computes the dot product of one vector with several others.
  !! @param[in] num_values_i the length of vector i.
  !! @param[in] indices_i the index value of the sparse vector i.
  !! @param[in] values_i the values of the sparse vector i.
  !! @param[in] num_values_j an array with the length of vectors j.
  !! @param[in] indices_j the indices of the vectors j.
  !! @param[in] values_j the values of the vectors j.
  !! @param[out] out_values the dot product values for each vector j.
  !! @param[in] comm the communicator to reduce along.
  SUBROUTINE DotAllHelper(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, out_values, comm)
    !! Parameters
    INTEGER, INTENT(IN) :: num_values_i
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    INTEGER, INTENT(INOUT) :: comm
    !! Local Variables
    INTEGER :: err
    INTEGER :: counter
    INTEGER :: inner_len_j

    !! Local Dot
    !$omp parallel private(inner_len_j)
    !$omp do
    DO counter = 1, SIZE(num_values_j)
       inner_len_j = num_values_j(counter)
       out_values(counter) = DotSparseVectors(indices_i(:num_values_i), &
            & values_i(:num_values_i), indices_j(:inner_len_j, counter), &
            & values_j(:inner_len_j, counter))
    END DO
    !$omp end do
    !$omp end parallel

    !! Reduce Over Processes
    CALL MPI_Allreduce(MPI_IN_PLACE, out_values, SIZE(num_values_j), &
         & MPINTREAL, MPI_SUM, comm, err)

  END SUBROUTINE DotAllHelper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DotAllPivoted(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, pivot_vector, num_local_pivots, out_values, comm)
    !! Parameters
    INTEGER, INTENT(IN) :: num_values_i
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    INTEGER, DIMENSION(:), INTENT(IN) :: pivot_vector
    INTEGER, INTENT(IN) :: num_local_pivots
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    INTEGER, INTENT(INOUT) :: comm
    !! Local Variables
    INTEGER :: err
    INTEGER :: counter
    INTEGER :: inner_len_j
    INTEGER :: local_pi_i

    !! Local Dot
    !$omp parallel private(inner_len_j, local_pi_i)
    !$omp do
    DO counter = 1, num_local_pivots
       local_pi_i = pivot_vector(counter)
       inner_len_j = num_values_j(local_pi_i)
       out_values(counter) = DotSparseVectors(indices_i(:num_values_i), &
            & values_i(:num_values_i), indices_j(:inner_len_j, local_pi_i), &
            & values_j(:inner_len_j, local_pi_i))
    END DO
    !$omp end do
    !$omp end parallel

    !! Reduce Over Processes
    CALL MPI_Allreduce(MPI_IN_PLACE, out_values, SIZE(num_values_j), &
         & MPINTREAL, MPI_SUM, comm, err)

  END SUBROUTINE DotAllPivoted
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to broadcast a sparse vector
  SUBROUTINE BroadcastVector(num_values, indices, values, root, comm)
    !! Parameters
    INTEGER, INTENT(INOUT) :: num_values
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    INTEGER, INTENT(IN) :: root
    INTEGER, INTENT(INOUT) :: comm
    !! Local
    INTEGER :: err

    CALL MPI_Bcast(num_values, 1, MPI_INT, root, comm, err)
    CALL MPI_Bcast(indices(:num_values), num_values, MPI_INT, root, comm, err)
    CALL MPI_Bcast(values(:num_values), num_values, MPINTREAL, root, comm, err)

  END SUBROUTINE BroadcastVector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to insert a value into a sparse vector.
  PURE SUBROUTINE AppendToVector(values_per, indices, values, insert_row, &
       & insert_value)
    !! Parameters
    INTEGER, INTENT(INOUT) :: values_per
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    INTEGER, INTENT(IN) :: insert_row
    REAL(NTREAL), INTENT(IN) :: insert_value

    values_per = values_per + 1
    indices(values_per) = insert_row
    values(values_per) = insert_value
  END SUBROUTINE AppendToVector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetPivot(AMat, start_index, pivot_vector, diag, index, &
       & local_pivots, num_local_pivots)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: AMat
    INTEGER, DIMENSION(:), INTENT(INOUT) :: pivot_vector
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: diag
    INTEGER, INTENT(IN) :: start_index
    INTEGER, INTENT(OUT) :: index
    INTEGER, DIMENSION(:), INTENT(INOUT) :: local_pivots
    INTEGER, INTENT(OUT) :: num_local_pivots
    !! Local Variables
    REAL(NTREAL) :: temp_diag
    DOUBLE PRECISION, DIMENSION(2) :: max_diag
    INTEGER :: pind
    INTEGER :: II
    INTEGER :: swap

    !! Search for the maximum pivot
    max_diag = [0, 0]
    DO II = start_index, AMat%actual_matrix_dimension
       pind = pivot_vector(II)
       IF (pind .GE. AMat%start_column .AND. pind .LT. AMat%end_column) THEN
          temp_diag = diag(pind - AMat%start_column + 1)
          IF (temp_diag .GT. max_diag(1)) THEN
             max_diag(1) = temp_diag
             max_diag(2) = II
          END IF
       END IF
    END DO
    CALL MPI_Allreduce(MPI_IN_PLACE, max_diag, 1, MPI_2DOUBLE_PRECISION, &
         & MPI_MAXLOC, row_comm, grid_error)
    index = INT(max_diag(2))

    swap = pivot_vector(index)
    pivot_vector(index) = pivot_vector(start_index)
    pivot_vector(start_index) = swap

    !! Determine local pivots
    num_local_pivots = 0
    DO II = start_index + 1, AMat%actual_matrix_dimension
       pind = pivot_vector(II)
       IF (pind .GE. AMat%start_column .AND. pind .LT. AMat%end_column) THEN
          num_local_pivots = num_local_pivots + 1
          local_pivots(num_local_pivots) = pind - AMat%start_column + 1
       END IF
    END DO
  END SUBROUTINE GetPivot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the vector holding the accumulated diagonal values
  SUBROUTINE ConstructDiag(AMat, dense_a, diag)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: AMat
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: dense_a
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: diag
    !! Local Variables
    INTEGER :: fill_counter
    INTEGER :: II, JJ, local_JJ, local_II
    INTEGER, DIMENSION(num_process_rows) :: diags_per_proc
    INTEGER, DIMENSION(num_process_rows) :: diag_displ

    diag = 0
    fill_counter = 0
    DO JJ = AMat%start_column, AMat%end_column - 1
       IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
          local_JJ = JJ - AMat%start_column + 1
          local_II = JJ - AMat%start_row + 1
          diag(local_JJ) = dense_a(local_II, local_JJ)
          fill_counter = fill_counter + 1
       END IF
    END DO
    diags_per_proc(my_row+1) = fill_counter
    !! Duplicate the diagonal entries along the process column (across rows)
    CALL MPI_Allgather(MPI_IN_PLACE, 1, MPI_INTEGER, diags_per_proc, 1, &
         & MPI_INTEGER, column_comm, grid_error)
    diag_displ(1) = 0
    DO II = 2, num_process_rows
       diag_displ(II) = diag_displ(II-1) + diags_per_proc(II-1)
    END DO
    CALL MPI_Allgatherv(MPI_IN_PLACE, diags_per_proc(my_row+1), MPINTREAL, &
         & diag, diags_per_proc, diag_displ, MPINTREAL, column_comm, grid_error)
  END SUBROUTINE ConstructDiag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ConstructRankLookup(AMat, col_root_lookup)
    !! Root Lookups
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: AMat
    INTEGER, DIMENSION(:), INTENT(INOUT) :: col_root_lookup
    !! Local Variables
    INTEGER, DIMENSION(num_process_columns) :: cols_per_proc
    INTEGER, DIMENSION(num_process_columns) :: d_cols_per_proc
    INTEGER :: II

    cols_per_proc(my_column+1) = AMat%local_columns
    CALL MPI_Allgather(MPI_IN_PLACE, 1, MPI_INTEGER, cols_per_proc, 1, &
         & MPI_INTEGER, row_comm, grid_error)
    d_cols_per_proc(1) = 0
    DO II = 2, num_process_columns
       d_cols_per_proc(II) = d_cols_per_proc(II-1) + cols_per_proc(II-1)
    END DO
    col_root_lookup(AMat%start_column:AMat%end_column - 1) = row_rank
    CALL MPI_Allgatherv(MPI_IN_PLACE, AMat%local_columns, MPI_INTEGER, &
         & col_root_lookup, cols_per_proc, d_cols_per_proc, MPI_INTEGER, &
         & row_comm, grid_error)

  END SUBROUTINE ConstructRankLookup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GatherMatrixColumn(local_matrix, column_matrix)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: local_matrix
    TYPE(SparseMatrix_t), INTENT(INOUT) :: column_matrix
    !! Local Variables
    TYPE(GatherHelper_t) :: gather_helper
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    TYPE(SparseMatrix_t) :: local_matrixT
    INTEGER :: mpi_err

    CALL TransposeSparseMatrix(local_matrix, local_matrixT)
    CALL GatherSizes(local_matrixT, column_comm, gather_helper)
    CALL MPI_Wait(gather_helper%size_request, mpi_status, mpi_err)
    CALL GatherAndComposeData(local_matrixT, column_comm, column_matrix, &
         & gather_helper)
    CALL MPI_Wait(gather_helper%outer_request,mpi_status,mpi_err)
    CALL MPI_Wait(gather_helper%inner_request,mpi_status,mpi_err)
    CALL MPI_Wait(gather_helper%data_request,mpi_status,mpi_err)
    CALL GatherAndComposeCleanup(local_matrixT, column_matrix, gather_helper)

    CALL DestructSparseMatrix(local_matrixT)
  END SUBROUTINE GatherMatrixColumn
END MODULE LinearSolversModule
