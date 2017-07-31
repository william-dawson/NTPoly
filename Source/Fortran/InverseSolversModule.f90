!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Inverse of a Matrix.
MODULE InverseSolversModule
  USE DataTypesModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ProcessGridModule
  USE TimerModule
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: Invert
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix.
  !! An implementation of Hotelling's method \cite palser1998canonical.
  !! @param[in] Mat1 Matrix 1.
  !! @param[out] InverseMat = Mat1^-1.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE Invert(Mat1, InverseMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: Mat1
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: InverseMat
    TYPE(IterativeSolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: TWO = 2.0
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Variables
    REAL(NTREAL) :: sigma
    TYPE(DistributedSparseMatrix_t) :: Temp1,Temp2,Identity
    TYPE(DistributedSparseMatrix_t) :: BalancedMat1
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL WriteCitation("palser1998canonical")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyDistributedSparseMatrix(InverseMat, Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Temp1, Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Temp2, Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Identity, Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(BalancedMat1, Mat1%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Mat1, BalancedMat1, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    ELSE
       CALL CopyDistributedSparseMatrix(Mat1,BalancedMat1)
    END IF
    CALL StopTimer("Load Balance")

    !! Compute Sigma
    CALL ComputeSigma(BalancedMat1,sigma)

    !! Create Inverse Guess
    CALL CopyDistributedSparseMatrix(BalancedMat1,InverseMat)
    CALL ScaleDistributedSparseMatrix(InverseMat,sigma)

    !! Iterate
    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    CALL StartTimer("I-Invert")
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       CALL DistributedGemm(InverseMat,BalancedMat1,Temp1, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

       !! Check if Converged
       CALL CopyDistributedSparseMatrix(Identity,Temp2)
       CALL IncrementDistributedSparseMatrix(Temp1,Temp2,REAL(-1.0,NTREAL))
       norm_value = DistributedSparseNorm(Temp2)

       CALL DestructDistributedSparseMatrix(Temp2)
       CALL DistributedGemm(Temp1,InverseMat,Temp2,alpha_in=REAL(-1.0,NTREAL), &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

       !! Save a copy of the last inverse matrix
       CALL CopyDistributedSparseMatrix(InverseMat,Temp1)

       CALL ScaleDistributedSparseMatrix(InverseMat,TWO)

       CALL IncrementDistributedSparseMatrix(Temp2,InverseMat, &
            & threshold_in=solver_parameters%threshold)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(InverseMat,min_size,max_size)
       sparsity = REAL(GetSize(InverseMat),KIND=NTREAL)/ &
            & (InverseMat%actual_matrix_dimension**2)
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
    END IF
    CALL StopTimer("I-Invert")

    !! Undo Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(InverseMat,InverseMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF
    CALL StopTimer("Load Balance")

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructDistributedSparseMatrix(Temp1)
    CALL DestructDistributedSparseMatrix(Temp2)
    CALL DestructDistributedSparseMatrix(BalancedMat1)
    CALL DestructDistributedMatrixMemoryPool(pool1)
  END SUBROUTINE Invert
END MODULE InverseSolversModule
