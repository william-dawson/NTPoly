!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Inverse of a Matrix.
MODULE InverseSolversModule
  USE DataTypesModule
  USE MatrixMemoryPoolPModule
  USE MatrixPSAlgebraModule
  USE MatrixPSModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: Invert
  PUBLIC :: PseudoInverse
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix.
  !! An implementation of Hotelling's method \cite palser1998canonical.
  !! @param[in] Mat1 Matrix 1.
  !! @param[out] InverseMat = Mat1^-1.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE Invert(Mat1, InverseMat, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: Mat1
    TYPE(Matrix_ps), INTENT(INOUT) :: InverseMat
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    REAL(NTREAL), PARAMETER :: TWO = 2.0
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Variables
    REAL(NTREAL) :: sigma
    TYPE(Matrix_ps) :: Temp1,Temp2,Identity
    TYPE(Matrix_ps) :: BalancedMat1
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool1

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL WriteCitation("palser1998canonical")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(InverseMat, Mat1)
    CALL ConstructEmptyMatrix(Temp1, Mat1)
    CALL ConstructEmptyMatrix(Temp2, Mat1)
    CALL ConstructEmptyMatrix(Identity, Mat1)
    CALL ConstructEmptyMatrix(BalancedMat1, Mat1)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Mat1, BalancedMat1, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    ELSE
       CALL CopyMatrix(Mat1,BalancedMat1)
    END IF
    CALL StopTimer("Load Balance")

    !! Compute Sigma
    CALL MatrixSigma(BalancedMat1,sigma)

    !! Create Inverse Guess
    CALL CopyMatrix(BalancedMat1,InverseMat)
    CALL ScaleMatrix(InverseMat,sigma)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
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

       CALL MatrixMultiply(InverseMat,BalancedMat1,Temp1, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

       !! Check if Converged
       CALL CopyMatrix(Identity,Temp2)
       CALL IncrementMatrix(Temp1,Temp2,REAL(-1.0,NTREAL))
       norm_value = MatrixNorm(Temp2)

       CALL DestructMatrix(Temp2)
       CALL MatrixMultiply(Temp1,InverseMat,Temp2,alpha_in=REAL(-1.0,NTREAL), &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

       !! Save a copy of the last inverse matrix
       CALL CopyMatrix(InverseMat,Temp1)

       CALL ScaleMatrix(InverseMat,TWO)

       CALL IncrementMatrix(Temp2,InverseMat, &
            & threshold_in=solver_parameters%threshold)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL PrintMatrixInformation(InverseMat)
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
    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(BalancedMat1)
    CALL DestructMatrixMemoryPool(pool1)
  END SUBROUTINE Invert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the pseudoinverse of a matrix.
  !! An implementation of Hotelling's method \cite palser1998canonical.
  !! We just change the convergence criteria.
  !! @param[in] Mat1 Matrix 1.
  !! @param[out] InverseMat = Mat1^-1.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE PseudoInverse(Mat1, InverseMat, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: Mat1
    TYPE(Matrix_ps), INTENT(INOUT) :: InverseMat
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    REAL(NTREAL), PARAMETER :: TWO = 2.0
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Variables
    REAL(NTREAL) :: sigma
    TYPE(Matrix_ps) :: Temp1,Temp2,Identity
    TYPE(Matrix_ps) :: BalancedMat1
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool1

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL WriteCitation("palser1998canonical")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(InverseMat, Mat1)
    CALL ConstructEmptyMatrix(Temp1, Mat1)
    CALL ConstructEmptyMatrix(Temp2, Mat1)
    CALL ConstructEmptyMatrix(Identity, Mat1)
    CALL ConstructEmptyMatrix(BalancedMat1, Mat1)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Mat1, BalancedMat1, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    ELSE
       CALL CopyMatrix(Mat1,BalancedMat1)
    END IF

    !! Compute Sigma
    CALL MatrixSigma(BalancedMat1,sigma)

    !! Create Inverse Guess
    CALL CopyMatrix(BalancedMat1,InverseMat)
    CALL ScaleMatrix(InverseMat,sigma)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       CALL MatrixMultiply(InverseMat,BalancedMat1,Temp1, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL MatrixMultiply(Temp1,InverseMat,Temp2,alpha_in=REAL(-1.0,NTREAL), &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

       !! Save a copy of the last inverse matrix
       CALL CopyMatrix(InverseMat,Temp1)

       CALL ScaleMatrix(InverseMat,TWO)
       CALL IncrementMatrix(Temp2,InverseMat, &
            & threshold_in=solver_parameters%threshold)

       !! Check if Converged
       CALL IncrementMatrix(InverseMat,Temp1,REAL(-1.0,NTREAL))
       norm_value = MatrixNorm(Temp1)

       !! Sometimes the first few values don't change so much, so that's why
       !! I added the outer counter check
       IF (norm_value .LE. solver_parameters%converge_diff .AND. &
            & outer_counter .GT. 3) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL PrintMatrixInformation(InverseMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(InverseMat,InverseMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(BalancedMat1)
    CALL DestructMatrixMemoryPool(pool1)
  END SUBROUTINE PseudoInverse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE InverseSolversModule
