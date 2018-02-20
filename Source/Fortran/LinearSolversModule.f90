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
    !! For Extracting To Global
    TYPE(TripletList_t) :: local_triplets
    TYPE(Triplet_t) :: temp
    !! Temporary Variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_index
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_values
    INTEGER :: recv_num_values
    INTEGER :: II, JJ, local_JJ, local_II, local_row
    REAL(NTREAL) :: dot_value, AJJ, insert_value, inverse_factor
    INTEGER :: root

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL MergeLocalBlocks(AMat, sparse_a)
    ALLOCATE(dense_a(sparse_a%rows, sparse_a%columns))
    dense_a = 0
    CALL ConstructDenseFromSparse(sparse_a, dense_a)

    !! Allocate space for L
    ALLOCATE(values_per_column_l(sparse_a%columns))
    ALLOCATE(index_l(sparse_a%rows, sparse_a%columns))
    ALLOCATE(values_l(sparse_a%rows, sparse_a%columns))
    values_per_column_l = 0

    !! Allocate space for a received column
    ALLOCATE(recv_index(sparse_a%rows))
    ALLOCATE(recv_values(sparse_a%rows))

    !! Main Loop
    DO JJ = 1, AMat%actual_matrix_dimension
       local_JJ = JJ - AMat%start_column + 1
       !! Dot Column JJ with Column JJ, Insert Value into L[J,J]
       inverse_factor = 0
       root = 0
       IF (JJ .GE. AMat%start_column .AND. JJ .LT. AMat%end_column) THEN
          dot_value = DotHelper(values_per_column_l(local_JJ), &
               & index_l(:,local_JJ), values_l(:,local_JJ), &
               & values_per_column_l(local_JJ), index_l(:,local_JJ), &
               & values_l(:,local_JJ), column_comm)
          IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
             local_row = JJ - AMat%start_row + 1
             AJJ = dense_a(local_row, local_JJ)
             insert_value = SQRT(AJJ - dot_value)
             inverse_factor = 1.0/insert_value
             !! Insert
             CALL AppendToVector(values_per_column_l(local_JJ), &
                  & index_l(:,local_JJ), values_l(:, local_JJ), local_row, &
                  & insert_value)
          END IF
          root = row_rank
       END IF
       !! Broadcast column JJ, and Inverse Factor
       CALL MPI_Allreduce(MPI_IN_PLACE, row_rank, 1, MPI_INTEGER, MPI_SUM, &
            & row_comm, grid_error)
       CALL BroadcastVector(recv_num_values, recv_index, recv_values, &
            & row_rank, row_comm)
       CALL MPI_Allreduce(MPI_IN_PLACE, inverse_factor, 1, MPINTREAL, MPI_SUM, &
            & within_slice_comm, grid_error)

       !! Loop over other columns
       DO II = JJ + 1, AMat%actual_matrix_dimension
          IF (II .GE. AMat%start_column .AND. II .LT. AMat%end_column) THEN
             local_II = II - AMat%start_column + 1
             dot_value = DotHelper(recv_num_values, recv_index, recv_values, &
                  values_per_column_l(local_II), index_l(:,local_II), &
                  & values_l(:,local_II), column_comm)
             IF (II .GE. AMat%start_row .AND. II .LT. AMat%end_row) THEN
                local_row = II - AMat%start_row + 1
                insert_value = inverse_factor * &
                     & (dense_a(local_row, local_II) - dot_value)
                CALL AppendToVector(values_per_column_l(local_II), &
                     & index_l(:,local_II), values_l(:, local_II), local_row, &
                     & insert_value)
             END IF
          END IF
       END DO
    END DO

    !! Extract The Result
    CALL ConstructTripletList(local_triplets)
    DO JJ = 1, SUM(values_per_column_l)
       !! note transpose
       temp%index_row = JJ + AMat%start_column - 1
       DO II = 1, values_per_column_l(JJ)
          !! note transpose
          temp%index_column = index_l(II,JJ) + AMat%start_row - 1
          temp%point_value = values_l(II,JJ)
          CALL AppendToTripletList(local_triplets, temp)
       END DO
    END DO

    CALL ConstructEmptyDistributedSparseMatrix(LMat, &
         & AMat%actual_matrix_dimension)
    CALL FillFromTripletList(LMat, local_triplets)

    !! Cleanup
    CALL DestructTripletList(local_triplets)
    DEALLOCATE(dense_a)
    DEALLOCATE(values_per_column_l)
    DEALLOCATE(index_l)
    DEALLOCATE(values_l)
    DEALLOCATE(recv_index)
    DEALLOCATE(recv_values)
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
    INTEGER, INTENT(IN), OPTIONAL :: rank_in
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! For Pivoting
    INTEGER, DIMENSION(:), ALLOCATABLE :: pivot_vector
    !! Temporary
    INTEGER :: i, j, counter

    !! Construct the pivot vector
    ALLOCATE(pivot_vector(AMat%actual_matrix_dimension))
    DO counter = 1, AMat%actual_matrix_dimension
       pivot_vector(counter) = counter
    END DO



    !! Cleanup
    DEALLOCATE(pivot_vector)
  END SUBROUTINE PivotedCholeskyDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine which computes sparse dot products across processors
  FUNCTION DotHelper(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, comm) RESULT(value)
    !! Parameters
    INTEGER, INTENT(IN) :: num_values_i, num_values_j
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i, indices_j
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i, values_j
    INTEGER, INTENT(INOUT) :: comm
    REAL(NTREAL) :: value
    !! Local Variables
    INTEGER :: err

    !! Local Dot
    value = DotSparseVectors(indices_i(:num_values_i), &
         & values_i(:num_values_i), indices_j(:num_values_j), &
         & values_j(:num_values_j))

    !! Reduce Over Processes
    CALL MPI_Allreduce(MPI_IN_PLACE, value, 1, MPINTREAL, MPI_SUM, comm, err)

  END FUNCTION DotHelper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to broadcast a sparse vector
  SUBROUTINE BroadcastVector(num_values, index, values, root, comm)
    !! Parameters
    INTEGER, INTENT(INOUT) :: num_values
    INTEGER, DIMENSION(:), INTENT(INOUT) :: index
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    INTEGER, INTENT(IN) :: root
    INTEGER, INTENT(INOUT) :: comm
    !! Local
    INTEGER :: rank
    INTEGER :: err

    CALL MPI_Bcast(num_values, 1, MPI_INT, root, comm, err)
    CALL MPI_Bcast(INDEX(:num_values), num_values, MPI_INT, root, comm, err)
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

    indices(values_per) = insert_row
    values(values_per) = insert_value
    values_per = values_per + 1
  END SUBROUTINE AppendToVector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine to append a value to a sparse matrix.
  SUBROUTINE AppendValue(matrix, column, row, value, insert_ptr)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matrix
    INTEGER, INTENT(IN) :: column
    INTEGER, INTENT(IN) :: row
    REAL(NTREAL), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: insert_ptr

    insert_ptr = insert_ptr + 1
    matrix%outer_index(column+1:) = matrix%outer_index(column+1) + 1
    matrix%inner_index(insert_ptr) = row
    matrix%values(insert_ptr) = value
  END SUBROUTINE AppendValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Symmetric Positive Definite matrix.
  !! This is a really naive implementation, that might be worth visiting.
  !! @param[in] AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] LMat the matrix computed.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE TempCholeskyDecomposition(AMat, LMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: AMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: LMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(TripletList_t) :: local_triplets
    TYPE(SparseMatrix_t) :: sparse_a, l_scratch, l_finished
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: dense_a
    !! Variables describing local matrix rows
    INTEGER :: local_row_j
    INTEGER :: local_row_i
    INTEGER :: local_column_j
    TYPE(SparseMatrix_t) :: row_j
    TYPE(SparseMatrix_t) :: row_i
    !! Temporary Variables
    REAL(NTREAL) :: local_dot
    REAL(NTREAL) :: global_dot
    REAL(NTREAL) :: insert_value
    REAL(NTREAL) :: inverse_factor
    INTEGER :: i, j, counter
    INTEGER :: insert_ptr
    INTEGER :: root

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL MergeLocalBlocks(AMat, sparse_a)
    ALLOCATE(dense_a(sparse_a%rows, sparse_a%columns))
    dense_a = 0
    CALL ConstructDenseFromSparse(sparse_a, dense_a)

    !! Make scratch space for the local matrix L
    CALL ConstructZeroSparseMatrix(l_scratch, sparse_a%rows, sparse_a%columns)
    DEALLOCATE(l_scratch%inner_index)
    DEALLOCATE(l_scratch%values)
    ALLOCATE(l_scratch%inner_index(l_scratch%rows*l_scratch%columns))
    l_scratch%inner_index = 0
    ALLOCATE(l_scratch%values(l_scratch%rows*l_scratch%columns))
    l_scratch%values = 0
    insert_ptr = 0

    !! Main Loop Over Rows
    DO j = 1, AMat%actual_matrix_dimension
       !! Diagonal Part
       local_row_j = j - AMat%start_row + 1
       local_column_j = j - AMat%start_column + 1
       IF (j .GE. AMat%start_row .AND. j .LT. AMat%end_row) THEN
          CALL ExtractRow(l_scratch, local_row_j, row_j)
          local_dot = DotSparseMatrix(row_j, row_j)
          CALL MPI_Allreduce(MPI_IN_PLACE, local_dot, 1, MPINTREAL, MPI_SUM, &
               & row_comm, grid_error)
          IF (j .GE. AMat%start_column .AND. j .LT. AMat%end_column) THEN
             local_dot = dense_a(local_row_j, local_column_j) - local_dot
             insert_value = SQRT(local_dot)
             inverse_factor = 1.0/insert_value
             CALL AppendValue(l_scratch, local_column_j, local_row_j, &
                  & insert_value, insert_ptr)
          END IF
          root = column_rank
       ELSE
          inverse_factor = 0
          root = 0
       END IF
       !! Broadcast Row J
       CALL MPI_Allreduce(MPI_IN_PLACE, root, 1, MPI_INTEGER, MPI_SUM, &
            & column_comm, grid_error)
       CALL BroadcastMatrix(row_j, column_comm, root)
       !! Broadcast Inverse Factor
       CALL MPI_Allreduce(MPI_IN_PLACE, inverse_factor, 1, MPINTREAL, MPI_SUM, &
            & column_comm, grid_error)
       !! Compute Dot Products
       DO i = j+1, AMat%actual_matrix_dimension
          IF (i .LT. AMat%start_row .OR. i .GE. AMat%end_row) CYCLE
          !! Extract Row I
          local_row_i = i - AMat%start_row + 1
          CALL ExtractRow(l_scratch, local_row_i, row_i)
          local_dot = DotSparseMatrix(row_i, row_j)
          CALL MPI_Allreduce(MPI_IN_PLACE, local_dot, 1, MPINTREAL, MPI_SUM, &
               & row_comm, grid_error)
          IF (j .GE. AMat%start_column .AND. j .LT. AMat%end_column) THEN
             local_dot = dense_a(local_row_i, local_column_j) - local_dot
             insert_value = inverse_factor*local_dot
             CALL AppendValue(l_scratch,local_column_j,local_row_i,insert_value,insert_ptr)
          END IF
       END DO
    END DO

    !! Trim extra memory off the scratch space
    CALL ConstructZeroSparseMatrix(l_finished, sparse_a%rows, sparse_a%columns)
    DEALLOCATE(l_finished%inner_index)
    DEALLOCATE(l_finished%values)
    ALLOCATE(l_finished%inner_index(insert_ptr))
    ALLOCATE(l_finished%values(insert_ptr))
    l_finished%outer_index = l_scratch%outer_index
    l_finished%inner_index = l_scratch%inner_index(:insert_ptr)
    l_finished%values = l_scratch%values(:insert_ptr)

    !! Finish by building the global L matrix
    IF (my_slice .EQ. 0) THEN
       CALL MatrixToTripletList(l_finished, local_triplets)
       DO counter = 1, local_triplets%CurrentSize
          local_triplets%data(counter)%index_row = &
               & local_triplets%data(counter)%index_row + AMat%start_row - 1
          local_triplets%data(counter)%index_column = &
               & local_triplets%data(counter)%index_column + AMat%start_column - 1
       END DO
    ELSE
       CALL ConstructTripletList(local_triplets)
    END IF
    CALL ConstructEmptyDistributedSparseMatrix(LMat, &
         & AMat%actual_matrix_dimension)
    CALL FillFromTripletList(LMat, local_triplets)

    !! Cleanup
    DEALLOCATE(dense_a)
    CALL DestructTripletList(local_triplets)
    CALL DestructSparseMatrix(sparse_a)
    CALL DestructSparseMatrix(l_scratch)
    CALL DestructSparseMatrix(l_finished)
  END SUBROUTINE TempCholeskyDecomposition
END MODULE LinearSolversModule
