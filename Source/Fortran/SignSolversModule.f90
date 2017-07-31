!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Matrix Sign Function.
MODULE SignSolversModule
  USE DataTypesModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixModule
  USE EigenBoundsModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ProcessGridModule
  USE TimerModule
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SignFunction
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function.
  !! @param[in] Mat1 the input matrix.
  !! @param[out] SignMat the sign of Mat1.
  !! @param[in] solver_parameters_in optional parameters for the routine.
  SUBROUTINE SignFunction(Mat1, SignMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: Mat1
    TYPE(DistributedSparseMatrix), INTENT(inout) :: SignMat
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: Identity
    TYPE(DistributedSparseMatrix) :: Temp1
    TYPE(DistributedSparseMatrix) :: Temp2
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Data
    REAL(NTREAL), PARAMETER :: alpha = 1.69770248526
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    REAL(NTREAL), PARAMETER :: THREE = 3.0
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: alpha_k
    REAL(NTREAL) :: xk
    REAL(NTREAL) :: norm_value
    INTEGER :: min_size, max_size
    INTEGER :: outer_counter
    REAL(NTREAL) :: sparsity

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Sign Function Solver")
       CALL EnterSubLog
       CALL WriteCitation("nicholas2008functions")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmpty(Identity, Mat1%actual_matrix_dimension)
    CALL ConstructEmpty(Temp1, Mat1%actual_matrix_dimension)
    CALL ConstructEmpty(Temp2, Mat1%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       !! Permute Matrices
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Mat1, SignMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    ELSE
       CALL CopyDistributedSparseMatrix(Mat1,SignMat)
    END IF
    CALL StopTimer("Load Balance")

    !! Initialize
    CALL GershgorinBounds(Mat1,e_min,e_max)
    xk = ABS(e_min/e_max)
    CALL ScaleDistributedSparseMatrix(SignMat,1.0/ABS(e_max))

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    iterate: DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       !! Update Scaling Factors
       alpha_k = MIN(SQRT(3.0/(1.0+xk+xk**2)), alpha)
       xk = 0.5*alpha_k*xk*(3.0-(alpha_k**2)*xk**2)

       CALL DistributedGemm(SignMat, SignMat, Temp1, &
            & alpha_in=-1.0*alpha_k**2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL IncrementDistributedSparseMatrix(Identity,Temp1,alpha_in=THREE)

       CALL DistributedGemm(SignMat, Temp1, Temp2, alpha_in=0.5*alpha_k, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

       CALL IncrementDistributedSparseMatrix(Temp2,SignMat,alpha_in=NEGATIVE_ONE)
       norm_value = DistributedSparseNorm(SignMat)
       CALL CopyDistributedSparseMatrix(Temp2,SignMat)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO iterate
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(SignMat,min_size,max_size)
       sparsity = GetSize(SignMat)/ &
            & (SignMat%actual_matrix_dimension**2)
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
    END IF

    !! Undo Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(SignMat,SignMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF
    CALL StopTimer("Load Balance")


    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructDistributedSparseMatrix(Temp1)
    CALL DestructDistributedSparseMatrix(Temp2)
    CALL DestructDistributedSparseMatrix(Identity)
  END SUBROUTINE SignFunction
END MODULE SignSolversModule
