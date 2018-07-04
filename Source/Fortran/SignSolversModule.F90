!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Matrix Sign Function.
MODULE SignSolversModule
  USE DataTypesModule
  USE MatrixMemoryPoolPModule
  USE MatrixPSAlgebraModule
  USE MatrixPSModule
  USE EigenBoundsModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SignFunction
  PUBLIC :: PolarDecomposition
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function.
  !! @param[in] Mat1 the input matrix.
  !! @param[out] SignMat the sign of Mat1.
  !! @param[in] solver_parameters_in optional parameters for the routine.
  SUBROUTINE SignFunction(Mat1, SignMat, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: Mat1
    TYPE(Matrix_ps), INTENT(INOUT) :: SignMat
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Sign Function Solver")
       CALL EnterSubLog
       CALL WriteCitation("nicholas2008functions")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    CALL CoreComputation(Mat1, SignMat, solver_parameters, .FALSE.)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE SignFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the polar decomposition of a matrix Mat1 = U*H.
  !! @param[in] Mat1 the input matrix.
  !! @param[out] Umat the unitary polar factor.
  !! @param[out] Hmat the hermitian matrix factor (optional).
  !! @param[in] solver_parameters_in optional parameters for the routine.
  SUBROUTINE PolarDecomposition(Mat1, Umat, Hmat, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: Mat1
    TYPE(Matrix_ps), INTENT(INOUT) :: Umat
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: Hmat
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    TYPE(Matrix_ps) :: UmatT

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Polar Decomposition Solver")
       CALL EnterSubLog
       CALL WriteCitation("nicholas2008functions")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    CALL CoreComputation(Mat1, Umat, solver_parameters, .TRUE.)

    IF (PRESENT(Hmat)) THEN
       CALL TransposeMatrix(Umat, UmatT)
       CALL MatrixMultiply(UmatT, Mat1, Hmat, &
            & threshold_in=solver_parameters%threshold)
       CALL DestructMatrix(UmatT)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE PolarDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This is the implementation routine for both the sign function and
  !! polar decomposition.
  SUBROUTINE CoreComputation(Mat1, OutMat, solver_parameters, needs_transpose)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: Mat1
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    TYPE(IterativeSolverParameters_t), INTENT(IN) :: solver_parameters
    LOGICAL, INTENT(IN) :: needs_transpose
    !! Local Matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: Temp1
    TYPE(Matrix_ps) :: Temp2
    TYPE(Matrix_ps) :: OutMatT
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Data
    REAL(NTREAL), PARAMETER :: alpha = 1.69770248526
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    REAL(NTREAL), PARAMETER :: THREE = 3.0
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: alpha_k
    REAL(NTREAL) :: xk
    REAL(NTREAL) :: norm_value
    INTEGER :: outer_counter

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Identity, Mat1)
    CALL ConstructEmptyMatrix(Temp1, Mat1)
    CALL ConstructEmptyMatrix(Temp2, Mat1)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       !! Permute Matrices
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Mat1, OutMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    ELSE
       CALL CopyMatrix(Mat1,OutMat)
    END IF
    CALL StopTimer("Load Balance")

    !! Initialize
    CALL GershgorinBounds(Mat1,e_min,e_max)
    xk = ABS(e_min/e_max)
    CALL ScaleMatrix(OutMat,1.0/ABS(e_max))

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

       IF (needs_transpose) THEN
          CALL TransposeMatrix(OutMat, OutMatT)
          CALL MatrixMultiply(OutMatT, OutMat, Temp1, &
               & alpha_in=-1.0*alpha_k**2, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       ELSE
          CALL MatrixMultiply(OutMat, OutMat, Temp1, &
               & alpha_in=-1.0*alpha_k**2, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       END IF
       CALL IncrementMatrix(Identity,Temp1,alpha_in=THREE)

       CALL MatrixMultiply(OutMat, Temp1, Temp2, alpha_in=0.5*alpha_k, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

       CALL IncrementMatrix(Temp2, OutMat, &
            & alpha_in=NEGATIVE_ONE)
       norm_value = MatrixNorm(OutMat)
       CALL CopyMatrix(Temp2,OutMat)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO iterate
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL PrintMatrixInformation(OutMat)
    END IF

    !! Undo Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat,OutMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF
    CALL StopTimer("Load Balance")

    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(OutMatT)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE CoreComputation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SignSolversModule
