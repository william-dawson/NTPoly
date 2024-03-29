!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Square Root of a Matrix.
MODULE SquareRootSolversModule
  USE ConvergenceMonitor, ONLY : ConstructMonitor, CheckConverged, AppendValue
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE EigenSolversModule, ONLY : DenseMatrixFunction
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteListElement, &
       & WriteHeader, WriteElement
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, &
       & IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, FillMatrixIdentity, PrintMatrixInformation
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: SquareRoot
  PUBLIC :: DenseSquareRoot
  PUBLIC :: InverseSquareRoot
  PUBLIC :: DenseInverseSquareRoot
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root of a matrix.
  SUBROUTINE SquareRoot(InputMat, OutputMat, solver_parameters_in, order_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The resulting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t),INTENT(IN),OPTIONAL :: solver_parameters_in
    !> Order of polynomial for calculation (default 5).
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    !! Call routine with the desired order.
    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, params, .FALSE.,&
            & order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, params, .FALSE.)
    END IF

    !! Cleanup
    CALL DestructSolverParameters(params)
  END SUBROUTINE SquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix square root function (dense version).
  SUBROUTINE DenseSquareRoot(Mat, OutputMat, solver_parameters_in)
    !> The matrix to compute the square root of.
    TYPE(Matrix_ps), INTENT(IN) :: Mat
    !> The computed matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Square Root Solver")
       CALL EnterSubLog
    END IF

    !! Apply
    CALL DenseMatrixFunction(Mat, OutputMat, SquareRootLambda, &
         & params)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructSolverParameters(params)
  END SUBROUTINE DenseSquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse square root of a matrix.
  SUBROUTINE InverseSquareRoot(InputMat, OutputMat, solver_parameters_in, &
       & order_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: InputMat
    !> The resulting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t),INTENT(IN),OPTIONAL :: solver_parameters_in
    !> Order of polynomial for calculation (default 5).
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    !! Call routine with the desired order.
    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, params, .TRUE., order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, params, .TRUE.)
    END IF

    !! Cleanup
    CALL DestructSolverParameters(params)

  END SUBROUTINE InverseSquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix inverse square root function (dense version).
  SUBROUTINE DenseInverseSquareRoot(Mat, OutputMat, solver_parameters_in)
    !> The matrix to compute the inverse square root of.
    TYPE(Matrix_ps), INTENT(IN)  :: Mat
    !> The computed matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Square Root Solver")
       CALL EnterSubLog
    END IF

    !! Apply
    CALL DenseMatrixFunction(Mat, OutputMat, InverseSquareRootLambda, params)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructSolverParameters(params)
  END SUBROUTINE DenseInverseSquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine picks the appropriate solver method
  SUBROUTINE SquareRootSelector(InputMat, OutputMat, params, &
       & compute_inverse, order_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The Matrix computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters about how to solve.
    TYPE(SolverParameters_t),INTENT(INOUT) :: params
    !> True if we are computing the inverse square root.
    LOGICAL, INTENT(IN) :: compute_inverse
    !> The polynomial degree to use (optional, default = 5)
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    INTEGER :: order

    IF (PRESENT(order_in)) THEN
       order = order_in
    ELSE
       order = 5
    END IF

    SELECT CASE(order)
    CASE(2)
       CALL NewtonSchultzISROrder2(InputMat, OutputMat, params, &
            & compute_inverse)
    CASE DEFAULT
       CALL NewtonSchultzISRTaylor(InputMat, OutputMat, params, &
            & order, compute_inverse)
    END SELECT

  END SUBROUTINE SquareRootSelector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root or inverse square root of a matrix.
  !> Based on the Newton-Schultz algorithm presented in: \cite jansik2007linear
  SUBROUTINE NewtonSchultzISROrder2(InMat, OutMat, params, compute_inverse)
    !> The matrix to compute
    TYPE(Matrix_ps), INTENT(IN)  :: InMat
    !> Mat^-1/2 or Mat^1/2.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(INOUT) :: params
    !> Whether to compute the inverse square root.
    LOGICAL, INTENT(IN) :: compute_inverse
    !! Local Variables
    REAL(NTREAL) :: lambda
    TYPE(Matrix_ps) :: X_k,T_k,Temp,Identity
    TYPE(Matrix_ps) :: SquareRootMat
    TYPE(Matrix_ps) :: InverseSquareRootMat
    !! Temporary Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: max_between
    INTEGER :: II
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: mpool

    !! Setup the monitor
    CALL ConstructMonitor(params%monitor, &
         & automatic_in = params%monitor_convergence, &
         & tight_cutoff_in=params%converge_diff)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("jansik2007linear")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(X_k, InMat)
    CALL ConstructEmptyMatrix(SquareRootMat, InMat)
    CALL ConstructEmptyMatrix(InverseSquareRootMat, InMat)
    CALL ConstructEmptyMatrix(T_k, InMat)
    CALL ConstructEmptyMatrix(Temp, InMat)
    CALL ConstructEmptyMatrix(Identity, InMat)
    CALL FillMatrixIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(InMat, e_min, e_max)
    max_between = MAX(ABS(e_min), ABS(e_max))
    lambda = 1.0 / max_between

    !! Initialize
    CALL FillMatrixIdentity(InverseSquareRootMat)
    CALL CopyMatrix(InMat,SquareRootMat)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat, SquareRootMat, &
            & params%BalancePermutation, memorypool_in = mpool)
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in = mpool)
       CALL PermuteMatrix(InverseSquareRootMat, InverseSquareRootMat, &
            & params%BalancePermutation, memorypool_in = mpool)
    END IF

    !! Iterate.
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    DO II = 1, params%max_iterations
       !! Compute X_k
       CALL MatrixMultiply(SquareRootMat, InverseSquareRootMat, X_k, &
            & threshold_in = params%threshold, memory_pool_in = mpool)
       CALL GershgorinBounds(X_k, e_min, e_max)
       max_between = MAX(ABS(e_min), ABS(e_max))
       lambda = 1.0 / max_between

       CALL ScaleMatrix(X_k, lambda)

       !! Check if Converged
       CALL CopyMatrix(Identity, Temp)
       CALL IncrementMatrix(X_k, Temp, &
            & alpha_in = -1.0_NTREAL)
       norm_value = MatrixNorm(Temp)

       !! Compute T_k
       CALL CopyMatrix(Identity, T_k)
       CALL ScaleMatrix(T_k, 3.0_NTREAL)
       CALL IncrementMatrix(X_k, T_k, &
            & alpha_in = -1.0_NTREAL)
       CALL ScaleMatrix(T_k, 0.5_NTREAL)

       !! Compute Z_k+1
       CALL CopyMatrix(InverseSquareRootMat, Temp)
       CALL MatrixMultiply(Temp, T_k, InverseSquareRootMat, &
            & threshold_in = params%threshold, memory_pool_in = mpool)
       CALL ScaleMatrix(InverseSquareRootMat, SQRT(lambda))

       !! Compute Y_k+1
       CALL CopyMatrix(SquareRootMat, Temp)
       CALL MatrixMultiply(T_k, Temp, SquareRootMat, &
            & threshold_in = params%threshold, memory_pool_in = mpool)
       CALL ScaleMatrix(SquareRootMat, SQRT(lambda))

       !! Check Exit Condition
       CALL AppendValue(params%monitor, norm_value)
       IF (CheckConverged(params%monitor, params%be_verbose)) EXIT

    END DO
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key = "Total Iterations", VALUE = II)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL CopyMatrix(InverseSquareRootMat, OutMat)
    ELSE
       CALL CopyMatrix(SquareRootMat, OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat, OutMat, &
            & params%BalancePermutation, memorypool_in = mpool)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructMatrix(Temp)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(SquareRootMat)
    CALL DestructMatrix(InverseSquareRootMat)
    CALL DestructMatrix(T_k)
    CALL DestructMatrixMemoryPool(mpool)
  END SUBROUTINE NewtonSchultzISROrder2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root or inverse square root of a matrix.
  !> Based on the Newton-Schultz algorithm with higher order polynomials.
  SUBROUTINE NewtonSchultzISRTaylor(InMat, OutMat, params, &
       & taylor_order, compute_inverse)
    !> Matrix to Compute
    TYPE(Matrix_ps), INTENT(IN) :: InMat
    !> Mat^-1/2 or Mat^1/2.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(INOUT) :: params
    !> Order of polynomial to use.
    INTEGER, INTENT(IN) :: taylor_order
    !> Whether to compute the inverse square root or not.
    LOGICAL, INTENT(IN) :: compute_inverse
    !! Local Variables
    REAL(NTREAL) :: lambda
    REAL(NTREAL) :: aa,bb,cc,dd
    REAL(NTREAL) :: a,b,c,d
    TYPE(Matrix_ps) :: X_k,Temp,Temp2,Identity
    TYPE(Matrix_ps) :: SquareRootMat
    TYPE(Matrix_ps) :: InverseSquareRootMat
    !! Temporary Variables
    REAL(NTREAL) :: e_min,e_max
    REAL(NTREAL) :: max_between
    INTEGER :: II
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: mpool

    !! Setup the monitor
    CALL ConstructMonitor(params%monitor, &
         & automatic_in = params%monitor_convergence, &
         & tight_cutoff_in=params%converge_diff)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("jansik2007linear")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(X_k, InMat)
    CALL ConstructEmptyMatrix(SquareRootMat, InMat)
    CALL ConstructEmptyMatrix(InverseSquareRootMat, InMat)
    CALL ConstructEmptyMatrix(Temp, InMat)
    IF (taylor_order == 5) THEN
       CALL ConstructEmptyMatrix(Temp2, InMat)
    END IF
    CALL ConstructEmptyMatrix(Identity, InMat)
    CALL FillMatrixIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(InMat, e_min, e_max)
    max_between = MAX(ABS(e_min), ABS(e_max))
    lambda = 1.0_NTREAL / max_between

    !! Initialize
    CALL FillMatrixIdentity(InverseSquareRootMat)
    CALL CopyMatrix(InMat, SquareRootMat)
    CALL ScaleMatrix(SquareRootMat, lambda)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat, SquareRootMat, &
            & params%BalancePermutation, memorypool_in = mpool)
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in = mpool)
       CALL PermuteMatrix(InverseSquareRootMat, InverseSquareRootMat, &
            & params%BalancePermutation, memorypool_in = mpool)
    END IF

    !! Iterate.
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    DO II = 1, params%max_iterations
       !! Compute X_k = Z_k * Y_k - I
       CALL MatrixMultiply(InverseSquareRootMat, SquareRootMat, X_k, &
            & threshold_in = params%threshold, memory_pool_in = mpool)
       CALL IncrementMatrix(Identity, X_k, -1.0_NTREAL)
       norm_value = MatrixNorm(X_k)

       SELECT CASE(taylor_order)
       CASE(3)
          !! Compute X_k^2
          CALL MatrixMultiply(X_k, X_k, Temp, &
               & threshold_in = params%threshold, memory_pool_in = mpool)

          !! X_k = I - 1/2 X_k + 3/8 X_k^2 + ...
          CALL ScaleMatrix(X_k, -0.5_NTREAL)
          CALL IncrementMatrix(Identity, X_k)
          CALL IncrementMatrix(Temp,X_k, 0.375_NTREAL)
       CASE(5)
          !! Compute p(x) = x^4 + A*x^3 + B*x^2 + C*x + D
          !! Scale to make coefficient of x^4 equal to 1
          aa = -40.0_NTREAL / 35.0_NTREAL
          bb = 48.0_NTREAL / 35.0_NTREAL
          cc = -64.0_NTREAL / 35.0_NTREAL
          dd = 128.0_NTREAL / 35.0_NTREAL

          !! The method of Knuth
          !! p = (z+x+b) * (z+c) + d
          !! z = x * (x+a)
          !! a = (A-1)/2
          !! b = B*(a+1) - C - a*(a+1)*(a+1)
          !! c = B - b - a*(a+1)
          !! d = D - b*c
          a = (aa - 1.0_NTREAL) / 2.0_NTREAL
          b = bb*(a + 1.0_NTREAL) - cc - a * (a + 1.0_NTREAL)**2
          c = bb - b - a * (a + 1.0_NTREAL)
          d = dd - b * c

          !! Compute Temp = z = x * (x+a)
          CALL MatrixMultiply(X_k, X_k, Temp, &
               & threshold_in = params%threshold, memory_pool_in = mpool)
          CALL IncrementMatrix(X_k, Temp, &
               & alpha_in = a)

          !! Compute Temp2 = z + x + b
          CALL CopyMatrix(Identity, Temp2)
          CALL ScaleMatrix(Temp2, b)
          CALL IncrementMatrix(X_k, Temp2)
          CALL IncrementMatrix(Temp, Temp2)

          !! Compute Temp = z + c
          CALL IncrementMatrix(Identity, Temp, c)

          !! Compute X_k = (z+x+b) * (z+c) + d = Temp2 * Temp + d
          CALL MatrixMultiply(Temp2, Temp, X_k, &
               & threshold_in = params%threshold, memory_pool_in = mpool)
          CALL IncrementMatrix(Identity, X_k, d)

          !! Scale back to the target coefficients
          CALL ScaleMatrix(X_k, 35.0_NTREAL / 128.0_NTREAL)
       END SELECT

       !! Compute Z_k+1 = Z_k * X_k
       CALL CopyMatrix(InverseSquareRootMat, Temp)
       CALL MatrixMultiply(X_k, Temp, InverseSquareRootMat, &
            & threshold_in = params%threshold, memory_pool_in = mpool)

       !! Compute Y_k+1 = X_k * Y_k
       CALL CopyMatrix(SquareRootMat, Temp)
       CALL MatrixMultiply(Temp, X_k, SquareRootMat, &
            & threshold_in = params%threshold, memory_pool_in = mpool)

       !! Check Exit Condition
       CALL AppendValue(params%monitor, norm_value)
       IF (CheckConverged(params%monitor, params%be_verbose)) EXIT

    END DO
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key = "Total Iterations", VALUE = II)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL ScaleMatrix(InverseSquareRootMat, SQRT(lambda))
       CALL CopyMatrix(InverseSquareRootMat, OutMat)
    ELSE
       CALL ScaleMatrix(SquareRootMat, 1.0_NTREAL / SQRT(lambda))
       CALL CopyMatrix(SquareRootMat, OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat, OutMat, &
            & params%BalancePermutation, memorypool_in = mpool)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructMatrix(X_k)
    CALL DestructMatrix(SquareRootMat)
    CALL DestructMatrix(InverseSquareRootMat)
    CALL DestructMatrix(Temp)
    IF (taylor_order == 5) THEN
       CALL DestructMatrix(Temp2)
    END IF
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(mpool)
  END SUBROUTINE NewtonSchultzISRTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prototypical square root function. 
  FUNCTION SquareRootLambda(val) RESULT(outval)
    REAL(KIND = NTREAL), INTENT(IN) :: val
    REAL(KIND = NTREAL) :: outval

    outval = SQRT(val)
  END FUNCTION SquareRootLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prototypical inverse square root function. 
  FUNCTION InverseSquareRootLambda(val) RESULT(outval)
    REAL(KIND = NTREAL), INTENT(IN) :: val
    REAL(KIND = NTREAL) :: outval

    outval = 1.0 / SQRT(val)
  END FUNCTION InverseSquareRootLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SquareRootSolversModule
