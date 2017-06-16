!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Solve the matrix equation AX = B
MODULE LinearSolversModule
  USE DataTypesModule
#ifdef NOBLOCK
  USE DistributedSparseMatrixModule
#else
  USE DistributedBlockedSparseMatrixModule
#endif
  USE DistributedMatrixMemoryPoolModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ProcessGridModule
  USE TimerModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: CGSolver
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Solve the matrix equation AX = B using conjugate gradient.
  !! @param[in] AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] XMat the solved for matrix X.
  !! @param[in] BMat the right hand side.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE CGSolver(AMat, XMat, BMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(IN)  :: AMat
    TYPE(DistributedSparseMatrix), INTENT(INOUT) :: XMat
    TYPE(DistributedSparseMatrix), INTENT(IN)  :: BMat
    TYPE(IterativeSolverParameters), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Variables
    TYPE(DistributedSparseMatrix) :: Identity
    TYPE(DistributedSparseMatrix) :: ABalanced
    TYPE(DistributedSparseMatrix) :: BBalanced
    TYPE(DistributedSparseMatrix) :: RMat, PMat, QMat
    TYPE(DistributedSparseMatrix) :: TempMat
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity
    REAL(NTREAL) :: top, bottom, new_top, step_size

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    !! Print out parameters
    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Linear Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="CG")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Setup all the matrices
    CALL ConstructEmpty(Identity, AMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)
    CALL ConstructEmpty(ABalanced, AMat%actual_matrix_dimension)
    CALL ConstructEmpty(BBalanced, AMat%actual_matrix_dimension)
    CALL ConstructEmpty(RMat, AMat%actual_matrix_dimension)
    CALL ConstructEmpty(PMat, AMat%actual_matrix_dimension)
    CALL ConstructEmpty(QMat, AMat%actual_matrix_dimension)
    CALL ConstructEmpty(TempMat, AMat%actual_matrix_dimension)

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
    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. IsRoot() .AND. &
            & outer_counter .GT. 1) THEN
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
       top = DotDistributedSparseMatrix(RMat,RMat)
       bottom = DotDistributedSparseMatrix(PMat,QMat)
       step_size = top/bottom

       !! Update
       CALL IncrementDistributedSparseMatrix(PMat, XMat, alpha_in=step_size)
       norm_value = step_size*DistributedSparseNorm(PMat)
       CALL IncrementDistributedSparseMatrix(QMat, RMat, alpha_in=-1.0*step_size)

       !! Update PMat
       new_top = DotDistributedSparseMatrix(RMat,RMat)
       step_size = new_top / top
       CALL ScaleDistributedSparseMatrix(PMat, step_size)
       CALL IncrementDistributedSparseMatrix(RMat, PMat)

    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(XMat,min_size,max_size)
       sparsity = GetSize(XMat)/ &
            & (XMat%actual_matrix_dimension**2)
       IF (IsRoot()) THEN
          CALL ExitSubLog
          CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
          CALL WriteHeader("Load_Balance")
          CALL EnterSubLog
          CALL WriteListElement(key="min_size", int_value_in=min_size)
          CALL WriteListElement(key="max_size", int_value_in=max_size)
          CALL ExitSubLog
          CALL WriteElement(key="Sparsity", float_value_in=sparsity)
       END IF
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
END MODULE LinearSolversModule