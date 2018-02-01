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
  USE ProcessGridModule
  USE SparseMatrixModule
  USE SparseVectorModule
  USE TimerModule
  USE TripletListModule
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
    TYPE(TripletList_t) :: local_triplets
    TYPE(SparseMatrix_t) :: sparse_a, sparse_l
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: dense_a
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_buffer
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: value_buffer
    !! Variables describing matrix row J
    INTEGER :: local_j
    INTEGER, DIMENSION(:), ALLOCATABLE :: index_j
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values_j
    !! Variables describing matrix row I
    INTEGER :: local_i
    INTEGER, DIMENSION(:), ALLOCATABLE :: index_i
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values_i
    !! Temporary Variables
    REAL(NTREAL) :: temp
    REAL(NTREAL) :: inverse_factor
    INTEGER :: i, j

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL GetTripletList(AMat, local_triplets)
    CALL ConstructFromTripletList(sparse_a, local_triplets, &
         & AMat%local_rows, AMat%local_columns)
    ALLOCATE(dense_a(sparse_a%rows,sparse_a%columns))
    dense_a = 0
    CALL ConstructDenseFromSparse(sparse_a, dense_a)

    !! Make scratch space for the local matrix L
    CALL ConstructZeroSparseMatrix(sparse_l, sparse_a%rows, sparse_a%columns)
    ALLOCATE(index_buffer(sparse_l%rows,sparse_l%columns))
    ALLOCATE(value_buffer(sparse_l%rows,sparse_l%columns))

    CALL PrintDistributedSparseMatrix(AMat)
    !! Main Loop Over Rows
    DO j = 1, AMat%actual_matrix_dimension
       !! Diagonal Part
       !! Extract row J
       local_j = j - AMat%start_row + 1
       
       !! Dot product, and then append diagonal value
       !!temp = DotSparseVectors(index_j, values_j, index_j, values_j)
       !! Scatter row j
       !! Compute Inverse Factor
       !! Compute Dot Products
       DO i = j+1, AMat%actual_matrix_dimension
         !! Extract Row I
         local_i = i - AMat%start_row + 1
         !! Dot product, adn then subtract from A
       END DO
    END DO

    !! Finish by building the global L matrix
    CALL MatrixToTripletList(sparse_l, local_triplets)
    CALL ConstructEmptyDistributedSparseMatrix(LMat, &
         & AMat%actual_matrix_dimension)
    CALL FillFromTripletList(LMat, local_triplets)

    !! Cleanup
    DEALLOCATE(dense_a)
    DEALLOCATE(index_buffer)
    DEALLOCATE(value_buffer)
    CALL DestructTripletList(local_triplets)
    CALL DestructSparseMatrix(sparse_a)
    CALL DestructSparseMatrix(sparse_l)
  END SUBROUTINE CholeskyDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Pivoted Cholesky Decomposition of a Symmetric Semi-Definite matrix.
  !! Pass it either a maximum rank and maximum error tolerance to check against
  !! @param[in] AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] LMat the matrix computed.
  !! @param[in] rank_in the target rank of the matrix.
  !! @param[out] rank_out the actual computed rank, if no target rank is
  !! specified
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE PivotedCholeskyDecomposition(AMat, LMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: AMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: LMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
  END SUBROUTINE PivotedCholeskyDecomposition
END MODULE LinearSolversModule
