!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Square Root of a Matrix.
MODULE SquareRootSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedMatrixMemoryPoolModule, ONLY : DistributedMatrixMemoryPool_t, &
       & DestructDistributedMatrixMemoryPool
  USE DistributedSparseMatrixAlgebraModule, ONLY : DistributedGemm, &
       & DistributedSparseNorm, IncrementDistributedSparseMatrix, &
       & ScaleDistributedSparseMatrix
  USE DistributedSparseMatrixModule, ONLY : DistributedSparseMatrix_t, &
       & ConstructEmptyDistributedSparseMatrix, CopyDistributedSparseMatrix, &
       & DestructDistributedSparseMatrix, FillDistributedIdentity, &
       & PrintMatrixInformation
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE IterativeSolversModule, ONLY : IterativeSolverParameters_t, &
       & PrintIterativeSolverParameters
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteHeader, WriteListElement, WriteCitation
  USE ProcessGridModule
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: SquareRoot
  PUBLIC :: InverseSquareRoot
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root of a matrix.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] solver_parameters_in parameters for the solver, optional.
  !! @param[in] order_in polynomial order for calculation (default 5).
  SUBROUTINE SquareRoot(InputMat, OutputMat, solver_parameters_in, order_in)
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: OutputMat
    TYPE(IterativeSolverParameters_t),INTENT(in),OPTIONAL :: &
         & solver_parameters_in
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    TYPE(IterativeSolverParameters_t) :: solver_parameters

    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .FALSE.,&
            & order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .FALSE.)
    END IF

  END SUBROUTINE SquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse square root of a matrix.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] solver_parameters_in parameters for the solver, optional.
  !! @param[in] order_in polynomial order for calculation (default 5).
  SUBROUTINE InverseSquareRoot(InputMat, OutputMat, solver_parameters_in, &
       & order_in)
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: OutputMat
    TYPE(IterativeSolverParameters_t),INTENT(in),OPTIONAL :: &
         & solver_parameters_in
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    TYPE(IterativeSolverParameters_t) :: solver_parameters

    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .TRUE.,&
            & order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .TRUE.)
    END IF

  END SUBROUTINE InverseSquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine picks the appropriate solver method
  !! @param[in] InputMat the matrix to compute.
  !! @param[inout] OutputMat the matrix computed.
  !! @param[in] solver_parameters parameters about how to solve.
  !! @param[in] compute_inverse true if we are computing the inverse square root
  !! @param[in] order_in the polynomial degree to use (optional, default=5)
  SUBROUTINE SquareRootSelector(InputMat, OutputMat, solver_parameters, &
       & compute_inverse, order_in)
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: OutputMat
    TYPE(IterativeSolverParameters_t),INTENT(in) :: solver_parameters
    LOGICAL, INTENT(IN) :: compute_inverse
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
  !! Based on the Newton-Schultz algorithm presented in: \cite jansik2007linear
  !! @param[in] Mat1 Matrix 1.
  !! @param[out] OutMat = Mat1^-1/2 or Mat1^1/2.
  !! @param[in] solver_parameters parameters for the solver
  !! @param[in] compute_inverse whether to compute the inverse square root.
  SUBROUTINE NewtonSchultzISROrder2(Mat1, OutMat, solver_parameters, &
       & compute_inverse)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: Mat1
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: OutMat
    TYPE(IterativeSolverParameters_t), INTENT(in) :: solver_parameters
    LOGICAL, INTENT(in) :: compute_inverse
    !! Local Variables
    REAL(NTREAL) :: lambda
    TYPE(DistributedSparseMatrix_t) :: X_k,T_k,Temp,Identity
    TYPE(DistributedSparseMatrix_t) :: SquareRootMat
    TYPE(DistributedSparseMatrix_t) :: InverseSquareRootMat
    !! Temporary Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: max_between
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(DistributedMatrixMemoryPool_t) :: pool1

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteCitation("jansik2007linear")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyDistributedSparseMatrix(X_k, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(SquareRootMat, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(InverseSquareRootMat, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(T_k, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Temp, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & Mat1%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(Mat1,e_min,e_max)
    max_between = MAX(ABS(e_min),ABS(e_max))
    lambda = 1.0/max_between

    !! Initialize
    CALL FillDistributedIdentity(InverseSquareRootMat)
    CALL CopyDistributedSparseMatrix(Mat1,SquareRootMat)

    !! Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat, SquareRootMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(InverseSquareRootMat, InverseSquareRootMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF
    CALL StopTimer("Load Balance")

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k
       CALL DistributedGemm(SquareRootMat,InverseSquareRootMat,X_k, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL GershgorinBounds(X_k,e_min,e_max)
       max_between = MAX(ABS(e_min),ABS(e_max))
       lambda = 1.0/max_between

       CALL ScaleDistributedSparseMatrix(X_k,lambda)

       !! Check if Converged
       CALL CopyDistributedSparseMatrix(Identity,Temp)
       CALL IncrementDistributedSparseMatrix(X_k,Temp,REAL(-1.0,NTREAL))
       norm_value = DistributedSparseNorm(Temp)

       !! Compute T_k
       CALL CopyDistributedSparseMatrix(Identity,T_k)
       CALL ScaleDistributedSparseMatrix(T_k,REAL(3.0,NTREAL))
       CALL IncrementDistributedSparseMatrix(X_k,T_k,REAL(-1.0,NTREAL))
       CALL ScaleDistributedSparseMatrix(T_k,REAL(0.5,NTREAL))

       !! Compute Z_k+1
       CALL CopyDistributedSparseMatrix(InverseSquareRootMat,Temp)
       CALL DistributedGemm(Temp,T_k,InverseSquareRootMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL ScaleDistributedSparseMatrix(InverseSquareRootMat,SQRT(lambda))

       !! Compute Y_k+1
       CALL CopyDistributedSparseMatrix(SquareRootMat, Temp)
       CALL DistributedGemm(T_k,Temp,SquareRootMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL ScaleDistributedSparseMatrix(SquareRootMat,SQRT(lambda))

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL CopyDistributedSparseMatrix(InverseSquareRootMat, OutMat)
    ELSE
       CALL CopyDistributedSparseMatrix(SquareRootMat, OutMat)
    END IF

    !! Undo Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat, OutMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF
    CALL StopTimer("Load Balance")

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructDistributedSparseMatrix(Temp)
    CALL DestructDistributedSparseMatrix(X_k)
    CALL DestructDistributedSparseMatrix(SquareRootMat)
    CALL DestructDistributedSparseMatrix(InverseSquareRootMat)
    CALL DestructDistributedSparseMatrix(T_k)
    CALL DestructDistributedMatrixMemoryPool(pool1)
  END SUBROUTINE NewtonSchultzISROrder2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root or inverse square root of a matrix.
  !! Based on the Newton-Schultz algorithm with higher order polynomials.
  !! @param[in] Mat1 Matrix 1.
  !! @param[out] OutMat = Mat1^-1/2 or Mat1^1/2.
  !! @param[in] solver_parameters parameters for the solver
  !! @param[in] compute_inverse whether to compute the inverse square root.
  SUBROUTINE NewtonSchultzISRTaylor(Mat1, OutMat, solver_parameters, &
       & taylor_order, compute_inverse)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: Mat1
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: OutMat
    TYPE(IterativeSolverParameters_t), INTENT(in) :: solver_parameters
    INTEGER, INTENT(in) :: taylor_order
    LOGICAL, INTENT(in) :: compute_inverse
    !! Local Variables
    REAL(NTREAL) :: lambda
    REAL(NTREAL) :: aa,bb,cc,dd
    REAL(NTREAL) :: a,b,c,d
    TYPE(DistributedSparseMatrix_t) :: X_k,Temp,Temp2,Identity
    TYPE(DistributedSparseMatrix_t) :: SquareRootMat
    TYPE(DistributedSparseMatrix_t) :: InverseSquareRootMat
    !! Temporary Variables
    REAL(NTREAL) :: e_min,e_max
    REAL(NTREAL) :: max_between
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(DistributedMatrixMemoryPool_t) :: pool1

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteCitation("jansik2007linear")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyDistributedSparseMatrix(X_k, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(SquareRootMat, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(InverseSquareRootMat, &
         & Mat1%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Temp, &
         & Mat1%actual_matrix_dimension)
    IF (taylor_order == 5) THEN
       CALL ConstructEmptyDistributedSparseMatrix(Temp2, &
            & Mat1%actual_matrix_dimension)
    END IF
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & Mat1%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(Mat1,e_min,e_max)
    max_between = MAX(ABS(e_min),ABS(e_max))
    lambda = 1.0_NTREAL/max_between

    !! Initialize
    CALL FillDistributedIdentity(InverseSquareRootMat)
    CALL CopyDistributedSparseMatrix(Mat1,SquareRootMat)
    CALL ScaleDistributedSparseMatrix(SquareRootMat,lambda)

    !! Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat,SquareRootMat, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
       CALL PermuteMatrix(Identity,Identity, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
       CALL PermuteMatrix(InverseSquareRootMat,InverseSquareRootMat, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
    END IF
    CALL StopTimer("Load Balance")

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k = Z_k * Y_k - I
       CALL DistributedGemm(InverseSquareRootMat,SquareRootMat,X_k, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)
       CALL IncrementDistributedSparseMatrix(Identity,X_k,-1.0_NTREAL)
       norm_value = DistributedSparseNorm(X_k)

       SELECT CASE(taylor_order)
       CASE(3)
          !! Compute X_k^2
          CALL DistributedGemm(X_k,X_k,Temp, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

          !! X_k = I - 1/2 X_k + 3/8 X_k^2 + ...
          CALL ScaleDistributedSparseMatrix(X_k,-0.5_NTREAL)
          CALL IncrementDistributedSparseMatrix(Identity,X_k)
          CALL IncrementDistributedSparseMatrix(Temp,X_k,0.375_NTREAL)
       CASE(5)
          !! Compute p(x) = x^4 + A*x^3 + B*x^2 + C*x + D
          !! Scale to make coefficient of x^4 equal to 1
          aa = -40.0_NTREAL/35.0_NTREAL
          bb = 48.0_NTREAL/35.0_NTREAL
          cc = -64.0_NTREAL/35.0_NTREAL
          dd = 128.0_NTREAL/35.0_NTREAL

          !! Knuth's method
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
          CALL DistributedGemm(X_k,X_k,Temp, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)
          CALL IncrementDistributedSparseMatrix(X_k,Temp,a)

          !! Compute Temp2 = z + x + b
          CALL CopyDistributedSparseMatrix(Identity,Temp2)
          CALL ScaleDistributedSparseMatrix(Temp2,b)
          CALL IncrementDistributedSparseMatrix(X_k,Temp2)
          CALL IncrementDistributedSparseMatrix(Temp,Temp2)

          !! Compute Temp = z + c
          CALL IncrementDistributedSparseMatrix(Identity,Temp,c)

          !! Compute X_k = (z+x+b) * (z+c) + d = Temp2 * Temp + d
          CALL DistributedGemm(Temp2,Temp,X_k, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)
          CALL IncrementDistributedSparseMatrix(Identity,X_k,d)

          !! Scale back to the target coefficients
          CALL ScaleDistributedSparseMatrix(X_k,35.0_NTREAL/128.0_NTREAL)
       END SELECT

       !! Compute Z_k+1 = Z_k * X_k
       CALL CopyDistributedSparseMatrix(InverseSquareRootMat,Temp)
       CALL DistributedGemm(X_k,Temp,InverseSquareRootMat, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

       !! Compute Y_k+1 = X_k * Y_k
       CALL CopyDistributedSparseMatrix(SquareRootMat,Temp)
       CALL DistributedGemm(Temp,X_k,SquareRootMat, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round",int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence",float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL ScaleDistributedSparseMatrix(InverseSquareRootMat,SQRT(lambda))
       CALL CopyDistributedSparseMatrix(InverseSquareRootMat,OutMat)
    ELSE
       CALL ScaleDistributedSparseMatrix(SquareRootMat,1.0_NTREAL/SQRT(lambda))
       CALL CopyDistributedSparseMatrix(SquareRootMat,OutMat)
    END IF

    !! Undo Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat,OutMat, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
    END IF
    CALL StopTimer("Load Balance")

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructDistributedSparseMatrix(X_k)
    CALL DestructDistributedSparseMatrix(SquareRootMat)
    CALL DestructDistributedSparseMatrix(InverseSquareRootMat)
    CALL DestructDistributedSparseMatrix(Temp)
    IF (taylor_order == 5) THEN
       CALL DestructDistributedSparseMatrix(Temp2)
    END IF
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedMatrixMemoryPool(pool1)
  END SUBROUTINE NewtonSchultzISRTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SquareRootSolversModule
