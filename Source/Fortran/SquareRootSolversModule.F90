!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Square Root of a Matrix.
MODULE SquareRootSolversModule
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
       & DestructSolverParameters
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
    TYPE(SolverParameters_t) :: solver_parameters

    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .FALSE.,&
            & order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .FALSE.)
    END IF

    !! Cleanup
    CALL DestructSolverParameters(solver_parameters)

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
    TYPE(SolverParameters_t) :: solver_parameters

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Square Root Solver")
       CALL EnterSubLog
    END IF

    !! Apply
    CALL DenseMatrixFunction(Mat, OutputMat, SquareRootLambda, &
         & solver_parameters)

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
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
    TYPE(SolverParameters_t) :: solver_parameters

    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .TRUE., &
            & order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .TRUE.)
    END IF

    !! Cleanup
    CALL DestructSolverParameters(solver_parameters)

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
    TYPE(SolverParameters_t) :: solver_parameters

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Square Root Solver")
       CALL EnterSubLog
    END IF

    !! Apply
    CALL DenseMatrixFunction(Mat, OutputMat, InverseSquareRootLambda, &
         & solver_parameters)

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE DenseInverseSquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine picks the appropriate solver method
  SUBROUTINE SquareRootSelector(InputMat, OutputMat, solver_parameters, &
       & compute_inverse, order_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The Matrix computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters about how to solve.
    TYPE(SolverParameters_t),INTENT(IN) :: solver_parameters
    !> True if we are computing the inverse square root.
    LOGICAL, INTENT(IN) :: compute_inverse
    !> The polynomial degree to use (optional, default=5)
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
       CALL NewtonSchultzISROrder2(InputMat, OutputMat, solver_parameters, &
            & compute_inverse)
    CASE DEFAULT
       CALL NewtonSchultzISRTaylor(InputMat, OutputMat, solver_parameters, &
            & order, compute_inverse)
    END SELECT

  END SUBROUTINE SquareRootSelector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root or inverse square root of a matrix.
  !> Based on the Newton-Schultz algorithm presented in: \cite jansik2007linear
  SUBROUTINE NewtonSchultzISROrder2(Mat, OutMat, solver_parameters, &
       & compute_inverse)
    !> The matrix to compute
    TYPE(Matrix_ps), INTENT(IN)  :: Mat
    !> Mat^-1/2 or Mat^1/2.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
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
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: mpool

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("jansik2007linear")
       CALL ExitSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(X_k, Mat)
    CALL ConstructEmptyMatrix(SquareRootMat, Mat)
    CALL ConstructEmptyMatrix(InverseSquareRootMat, Mat)
    CALL ConstructEmptyMatrix(T_k, Mat)
    CALL ConstructEmptyMatrix(Temp, Mat)
    CALL ConstructEmptyMatrix(Identity, Mat)
    CALL FillMatrixIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(Mat,e_min,e_max)
    max_between = MAX(ABS(e_min),ABS(e_max))
    lambda = 1.0/max_between

    !! Initialize
    CALL FillMatrixIdentity(InverseSquareRootMat)
    CALL CopyMatrix(Mat,SquareRootMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat, SquareRootMat, &
            & solver_parameters%BalancePermutation, memorypool_in=mpool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=mpool)
       CALL PermuteMatrix(InverseSquareRootMat, InverseSquareRootMat, &
            & solver_parameters%BalancePermutation, memorypool_in=mpool)
    END IF

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k
       CALL MatrixMultiply(SquareRootMat,InverseSquareRootMat,X_k, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=mpool)
       CALL GershgorinBounds(X_k,e_min,e_max)
       max_between = MAX(ABS(e_min),ABS(e_max))
       lambda = 1.0/max_between

       CALL ScaleMatrix(X_k,lambda)

       !! Check if Converged
       CALL CopyMatrix(Identity,Temp)
       CALL IncrementMatrix(X_k,Temp,REAL(-1.0,NTREAL))
       norm_value = MatrixNorm(Temp)

       !! Compute T_k
       CALL CopyMatrix(Identity,T_k)
       CALL ScaleMatrix(T_k,REAL(3.0,NTREAL))
       CALL IncrementMatrix(X_k,T_k,REAL(-1.0,NTREAL))
       CALL ScaleMatrix(T_k,REAL(0.5,NTREAL))

       !! Compute Z_k+1
       CALL CopyMatrix(InverseSquareRootMat,Temp)
       CALL MatrixMultiply(Temp,T_k,InverseSquareRootMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=mpool)
       CALL ScaleMatrix(InverseSquareRootMat,SQRT(lambda))

       !! Compute Y_k+1
       CALL CopyMatrix(SquareRootMat, Temp)
       CALL MatrixMultiply(T_k,Temp,SquareRootMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=mpool)
       CALL ScaleMatrix(SquareRootMat,SQRT(lambda))

       IF (solver_parameters%be_verbose) THEN
          CALL WriteElement(key="Convergence", VALUE=norm_value)
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", VALUE=outer_counter)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL CopyMatrix(InverseSquareRootMat, OutMat)
    ELSE
       CALL CopyMatrix(SquareRootMat, OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat, OutMat, &
            & solver_parameters%BalancePermutation, memorypool_in=mpool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
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
  SUBROUTINE NewtonSchultzISRTaylor(Mat, OutMat, solver_parameters, &
       & taylor_order, compute_inverse)
    !> Matrix to Compute
    TYPE(Matrix_ps), INTENT(IN)  :: Mat
    !> Mat^-1/2 or Mat^1/2.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
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
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: mpool

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("jansik2007linear")
       CALL ExitSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(X_k, Mat)
    CALL ConstructEmptyMatrix(SquareRootMat, Mat)
    CALL ConstructEmptyMatrix(InverseSquareRootMat, Mat)
    CALL ConstructEmptyMatrix(Temp, Mat)
    IF (taylor_order == 5) THEN
       CALL ConstructEmptyMatrix(Temp2, Mat)
    END IF
    CALL ConstructEmptyMatrix(Identity, Mat)
    CALL FillMatrixIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(Mat,e_min,e_max)
    max_between = MAX(ABS(e_min),ABS(e_max))
    lambda = 1.0_NTREAL/max_between

    !! Initialize
    CALL FillMatrixIdentity(InverseSquareRootMat)
    CALL CopyMatrix(Mat,SquareRootMat)
    CALL ScaleMatrix(SquareRootMat,lambda)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat,SquareRootMat, &
            & solver_parameters%BalancePermutation,memorypool_in=mpool)
       CALL PermuteMatrix(Identity,Identity, &
            & solver_parameters%BalancePermutation,memorypool_in=mpool)
       CALL PermuteMatrix(InverseSquareRootMat,InverseSquareRootMat, &
            & solver_parameters%BalancePermutation,memorypool_in=mpool)
    END IF

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k = Z_k * Y_k - I
       CALL MatrixMultiply(InverseSquareRootMat,SquareRootMat,X_k, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=mpool)
       CALL IncrementMatrix(Identity,X_k,-1.0_NTREAL)
       norm_value = MatrixNorm(X_k)

       SELECT CASE(taylor_order)
       CASE(3)
          !! Compute X_k^2
          CALL MatrixMultiply(X_k,X_k,Temp, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=mpool)

          !! X_k = I - 1/2 X_k + 3/8 X_k^2 + ...
          CALL ScaleMatrix(X_k,-0.5_NTREAL)
          CALL IncrementMatrix(Identity,X_k)
          CALL IncrementMatrix(Temp,X_k,0.375_NTREAL)
       CASE(5)
          !! Compute p(x) = x^4 + A*x^3 + B*x^2 + C*x + D
          !! Scale to make coefficient of x^4 equal to 1
          aa = -40.0_NTREAL/35.0_NTREAL
          bb = 48.0_NTREAL/35.0_NTREAL
          cc = -64.0_NTREAL/35.0_NTREAL
          dd = 128.0_NTREAL/35.0_NTREAL

          !! The method of Knuth
          !! p = (z+x+b) * (z+c) + d
          !! z = x * (x+a)
          !! a = (A-1)/2
          !! b = B*(a+1) - C - a*(a+1)*(a+1)
          !! c = B - b - a*(a+1)
          !! d = D - b*c
          a = (aa-1.0_NTREAL)/2.0_NTREAL
          b = bb*(a+1.0_NTREAL)-cc-a*(a+1.0_NTREAL)**2
          c = bb-b-a*(a+1.0_NTREAL)
          d = dd-b*c

          !! Compute Temp = z = x * (x+a)
          CALL MatrixMultiply(X_k,X_k,Temp, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=mpool)
          CALL IncrementMatrix(X_k,Temp,a)

          !! Compute Temp2 = z + x + b
          CALL CopyMatrix(Identity,Temp2)
          CALL ScaleMatrix(Temp2,b)
          CALL IncrementMatrix(X_k,Temp2)
          CALL IncrementMatrix(Temp,Temp2)

          !! Compute Temp = z + c
          CALL IncrementMatrix(Identity,Temp,c)

          !! Compute X_k = (z+x+b) * (z+c) + d = Temp2 * Temp + d
          CALL MatrixMultiply(Temp2,Temp,X_k, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=mpool)
          CALL IncrementMatrix(Identity,X_k,d)

          !! Scale back to the target coefficients
          CALL ScaleMatrix(X_k,35.0_NTREAL/128.0_NTREAL)
       END SELECT

       !! Compute Z_k+1 = Z_k * X_k
       CALL CopyMatrix(InverseSquareRootMat,Temp)
       CALL MatrixMultiply(X_k,Temp,InverseSquareRootMat, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=mpool)

       !! Compute Y_k+1 = X_k * Y_k
       CALL CopyMatrix(SquareRootMat,Temp)
       CALL MatrixMultiply(Temp,X_k,SquareRootMat, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=mpool)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", VALUE=outer_counter)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL ScaleMatrix(InverseSquareRootMat,SQRT(lambda))
       CALL CopyMatrix(InverseSquareRootMat,OutMat)
    ELSE
       CALL ScaleMatrix(SquareRootMat,1.0_NTREAL/SQRT(lambda))
       CALL CopyMatrix(SquareRootMat,OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat,OutMat, &
            & solver_parameters%BalancePermutation,memorypool_in=mpool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
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
  !> Prototypical sine function. 
  SUBROUTINE SquareRootLambda(index, val)
    !> The index of the eigenvalue
    INTEGER, INTENT(IN) :: index
    !> The actual value of an element.
    REAL(KIND=NTREAL), INTENT(INOUT) :: val

    val = SQRT(val)
  END SUBROUTINE SquareRootLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prototypical sine function. 
  SUBROUTINE InverseSquareRootLambda(index, val)
    !> The index of the eigenvalue
    INTEGER, INTENT(IN) :: index
    !> The actual value of an element.
    REAL(KIND=NTREAL), INTENT(INOUT) :: val

    val = 1.0/SQRT(val)
  END SUBROUTINE InverseSquareRootLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SquareRootSolversModule
