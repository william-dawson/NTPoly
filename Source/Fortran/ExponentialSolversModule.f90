!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Matrix Exponentials and Logarithms.
MODULE ExponentialSolversModule
  USE ChebyshevSolversModule, ONLY : ChebyshevPolynomial_t, &
       & ChebyshevCompute, ConstructChebyshevPolynomial, &
       & DestructChebyshevPolynomial, FactorizedChebyshevCompute, &
       & SetChebyshevCoefficient
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixAlgebraModule, ONLY : DistributedGemm, &
       & DistributedSparseNorm, IncrementDistributedSparseMatrix, &
       & ScaleDistributedSparseMatrix
  USE DistributedMatrixMemoryPoolModule, ONLY : DistributedMatrixMemoryPool_t
  USE DistributedSparseMatrixModule, ONLY : DistributedSparseMatrix_t, &
       & ConstructEmptyDistributedSparseMatrix, CopyDistributedSparseMatrix, &
       & DestructDistributedSparseMatrix, FillDistributedIdentity, &
       & PrintMatrixInformation
  USE EigenBoundsModule, ONLY : GershgorinBounds, PowerBounds
  USE FixedSolversModule, ONLY : FixedSolverParameters_t, &
       & PrintFixedSolverParameters
  USE IterativeSolversModule, ONLY : IterativeSolverParameters_t, &
       & PrintIterativeSolverParameters
  USE LinearSolversModule, ONLY : CGSolver
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteHeader, WriteElement, WriteListElement, &
       & EnterSubLog, ExitSubLog
  USE ParameterConverterModule, ONLY : ConvertFixedToIterative
  USE ProcessGridModule
  USE RootSolversModule, ONLY : ComputeRoot
  USE SquareRootSolversModule, ONLY : SquareRoot
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: ComputeExponential
  PUBLIC :: ComputeExponentialPade
  PUBLIC :: ComputeExponentialTaylor
  PUBLIC :: ComputeLogarithm
  PUBLIC :: ComputeLogarithmTaylor
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = exp(InputMat)
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeExponential(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: OutputMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    TYPE(FixedSolverParameters_t) :: sub_solver_parameters
    TYPE(IterativeSolverParameters_t) :: i_sub_solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: ScaledMat
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! For Chebyshev Expansion
    TYPE(ChebyshevPolynomial_t) :: polynomial
    !! Local Variables
    REAL(NTREAL) :: spectral_radius
    REAL(NTREAL) :: sigma_val
    INTEGER :: sigma_counter
    INTEGER :: counter

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF
    sub_solver_parameters = solver_parameters
    CALL ConvertFixedToIterative(solver_parameters, i_sub_solver_parameters)
    i_sub_solver_parameters%max_iterations = 10

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Chebyshev")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    CALL ConstructEmptyDistributedSparseMatrix(OutputMat, &
         & InputMat%actual_matrix_dimension)

    !! Scale the matrix
    CALL PowerBounds(InputMat,spectral_radius,i_sub_solver_parameters)
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 1.0)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO
    CALL CopyDistributedSparseMatrix(InputMat, ScaledMat)
    CALL ScaleDistributedSparseMatrix(ScaledMat,1.0/sigma_val)
    sub_solver_parameters%threshold = sub_solver_parameters%threshold/sigma_val

    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Sigma", float_value_in=sigma_val)
    END IF

    !! Expand Chebyshev Series
    CALL ConstructChebyshevPolynomial(polynomial,16)
    CALL SetChebyshevCoefficient(polynomial,1,1.266065877752007e+00_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,2,1.130318207984970e+00_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,3,2.714953395340771e-01_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,4,4.433684984866504e-02_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,5,5.474240442092110e-03_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,6,5.429263119148932e-04_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,7,4.497732295351912e-05_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,8,3.198436462630565e-06_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,9,1.992124801999838e-07_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,10,1.103677287249654e-08_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,11,5.505891628277851e-10_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,12,2.498021534339559e-11_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,13,1.038827668772902e-12_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,14,4.032447357431817e-14_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,15,2.127980007794583e-15_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,16,-1.629151584468762e-16_NTREAL)

    CALL ChebyshevCompute(ScaledMat,OutputMat,polynomial,sub_solver_parameters)
    !CALL FactorizedChebyshevCompute(ScaledMat,OutputMat,polynomial, &
    !     & sub_solver_parameters)

    !! Undo the scaling by squaring at the end.
    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    DO counter=1,sigma_counter-1
       CALL DistributedGemm(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparseMatrix(TempMat,OutputMat)
    END DO

    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(OutputMat)
    END IF

    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructChebyshevPolynomial(polynomial)
    CALL DestructDistributedSparseMatrix(ScaledMat)
    CALL DestructDistributedSparseMatrix(TempMat)
  END SUBROUTINE ComputeExponential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix using a pade approximation.
  !! Be warned, the pade method can result in a lot of intermediate fill.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = exp(InputMat)
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeExponentialPade(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: OutputMat
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Solver Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    TYPE(IterativeSolverParameters_t) :: sub_solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: ScaledMat
    TYPE(DistributedSparseMatrix_t) :: IdentityMat
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedSparseMatrix_t) :: B1, B2, B3
    TYPE(DistributedSparseMatrix_t) :: P1, P2
    TYPE(DistributedSparseMatrix_t) :: LeftMat, RightMat
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Variables
    REAL(NTREAL) :: spectral_radius
    REAL(NTREAL) :: sigma_val
    INTEGER :: sigma_counter
    INTEGER :: counter

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Pade")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Setup
    CALL ConstructEmptyDistributedSparseMatrix(IdentityMat, &
         & InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(IdentityMat)

    !! Scale the matrix
    spectral_radius = DistributedSparseNorm(InputMat)
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 1.0)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO
    CALL CopyDistributedSparseMatrix(InputMat, ScaledMat)
    CALL ScaleDistributedSparseMatrix(ScaledMat,1.0/sigma_val)
    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Sigma", float_value_in=sigma_val)
       CALL WriteElement(key="Scaling_Steps", int_value_in=sigma_counter)
    END IF

    !! Sub Solver Parameters
    sub_solver_parameters = solver_parameters
    sub_solver_parameters%threshold = sub_solver_parameters%threshold/sigma_val

    !! Power Matrices
    CALL DistributedGemm(ScaledMat, ScaledMat, B1, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)
    CALL DistributedGemm(B1, B1, B2, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)
    CALL DistributedGemm(B2, B2, B3, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)

    !! Polynomials - 1
    CALL CopyDistributedSparseMatrix(IdentityMat, P1)
    CALL ScaleDistributedSparseMatrix(P1,17297280.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(B1, P1, 1995840.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(B2, P1, 25200.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(B3, P1, 56.0_NTREAL)
    !! Polynomials - 2
    CALL CopyDistributedSparseMatrix(IdentityMat, TempMat)
    CALL ScaleDistributedSparseMatrix(TempMat,8648640.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(B1, TempMat, 277200.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(B2, TempMat, 1512.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(B3, TempMat)
    CALL DistributedGemm(ScaledMat, TempMat, P2, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)

    !! Left and Right
    CALL CopyDistributedSparseMatrix(P1, LeftMat)
    CALL IncrementDistributedSparseMatrix(P2, LeftMat, -1.0_NTREAL)
    CALL CopyDistributedSparseMatrix(P1, RightMat)
    CALL IncrementDistributedSparseMatrix(P2, RightMat, 1.0_NTREAL)

    CALL CGSolver(LeftMat, OutputMat, RightMat, sub_solver_parameters)

    !! Undo the scaling by squaring at the end.
    DO counter=1,sigma_counter-1
       CALL DistributedGemm(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparseMatrix(TempMat,OutputMat)
    END DO

    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructDistributedSparseMatrix(ScaledMat)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(B1)
    CALL DestructDistributedSparseMatrix(B2)
    CALL DestructDistributedSparseMatrix(B3)
    CALL DestructDistributedSparseMatrix(P1)
    CALL DestructDistributedSparseMatrix(P2)
    CALL DestructDistributedSparseMatrix(LeftMat)
    CALL DestructDistributedSparseMatrix(RightMat)
  END SUBROUTINE ComputeExponentialPade
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix using a taylor series expansion.
  !! This is only really useful if you have a very small spectrum, because
  !! quite a bit of scaling is required.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = exp(InputMat)
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeExponentialTaylor(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: OutputMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    TYPE(IterativeSolverParameters_t) :: i_sub_solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: ScaledMat
    TYPE(DistributedSparseMatrix_t) :: Ak
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Variables
    REAL(NTREAL) :: spectral_radius
    REAL(NTREAL) :: sigma_val
    REAL(NTREAL) :: taylor_denom
    INTEGER :: sigma_counter
    INTEGER :: counter

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF
    CALL ConvertFixedToIterative(solver_parameters, i_sub_solver_parameters)
    i_sub_solver_parameters%max_iterations = 10

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Taylor")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL PowerBounds(InputMat,spectral_radius,i_sub_solver_parameters)

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 3.0e-8)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL CopyDistributedSparseMatrix(InputMat, ScaledMat)
    CALL ScaleDistributedSparseMatrix(ScaledMat,1.0/sigma_val)

    CALL ConstructEmptyDistributedSparseMatrix(OutputMat, &
         & InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(OutputMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(ScaledMat, ScaledMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Expand Taylor Series
    taylor_denom = 1.0
    CALL CopyDistributedSparseMatrix(OutputMat, Ak)
    DO counter=1,10
       taylor_denom = taylor_denom * counter
       CALL DistributedGemm(Ak,ScaledMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparseMatrix(TempMat,Ak)
       CALL IncrementDistributedSparseMatrix(Ak,OutputMat)
    END DO

    DO counter=1,sigma_counter-1
       CALL DistributedGemm(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparseMatrix(TempMat,OutputMat)
    END DO

    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructDistributedSparseMatrix(ScaledMat)
    CALL DestructDistributedSparseMatrix(Ak)
    CALL DestructDistributedSparseMatrix(TempMat)
  END SUBROUTINE ComputeExponentialTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the logarithm of a matrix.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = log(InputMat)
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeLogarithm(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: OutputMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: ScaledMat
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedSparseMatrix_t) :: IdentityMat
    !! For Chebyshev Expansion
    TYPE(ChebyshevPolynomial_t) :: polynomial
    !! Local Variables
    TYPE(IterativeSolverParameters_t) :: i_sub_solver_parameters
    TYPE(IterativeSolverParameters_t) :: p_sub_solver_parameters
    TYPE(FixedSolverParameters_t) :: f_sub_solver_parameters
    REAL(NTREAL) :: spectral_radius
    INTEGER :: sigma_val
    INTEGER :: sigma_counter

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF
    CALL ConvertFixedToIterative(solver_parameters, i_sub_solver_parameters)
    CALL ConvertFixedToIterative(solver_parameters, p_sub_solver_parameters)
    p_sub_solver_parameters%max_iterations=16
    f_sub_solver_parameters = solver_parameters

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Logarithm Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Chebyshev")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Setup
    CALL ConstructEmptyDistributedSparseMatrix(IdentityMat, &
         & InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(IdentityMat)

    !! Copy to a temporary matrix for scaling.
    CALL CopyDistributedSparseMatrix(InputMat,ScaledMat)

    !! Compute The Scaling Factor
    sigma_val = 1
    sigma_counter = 1
    CALL PowerBounds(InputMat,spectral_radius,p_sub_solver_parameters)
    DO WHILE (spectral_radius .GT. SQRT(2.0))
       spectral_radius = SQRT(spectral_radius)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Sigma", int_value_in=sigma_val)
    END IF
    f_sub_solver_parameters%threshold = &
         & f_sub_solver_parameters%threshold/REAL(2**(sigma_counter-1),NTREAL)
    CALL ComputeRoot(InputMat, ScaledMat, sigma_val, i_sub_solver_parameters)

    !! Shift Scaled Matrix
    CALL IncrementDistributedSparseMatrix(IdentityMat,ScaledMat, &
         & alpha_in=REAL(-1.0,NTREAL))

    !! Expand Chebyshev Series
    CALL ConstructChebyshevPolynomial(polynomial,32)
    CALL SetChebyshevCoefficient(polynomial,1,-0.485101351704_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,2,1.58828112379_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,3,-0.600947731795_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,4,0.287304748177_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,5,-0.145496447103_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,6,0.0734013668818_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,7,-0.0356277942958_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,8,0.0161605505166_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,9,-0.0066133591188_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,10,0.00229833505456_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,11,-0.000577804103964_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,12,2.2849332964e-05_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,13,8.37426826403e-05_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,14,-6.10822859027e-05_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,15,2.58132364523e-05_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,16,-5.87577322647e-06_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,17,-8.56711062722e-07_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,18,1.52066488969e-06_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,19,-7.12760496253e-07_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,20,1.23102245249e-07_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,21,6.03168259043e-08_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,22,-5.1865499826e-08_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,23,1.43185107512e-08_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,24,2.58449717089e-09_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,25,-3.73189861771e-09_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,26,1.18469334815e-09_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,27,1.51569931066e-10_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,28,-2.89595999673e-10_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,29,1.26720668874e-10_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,30,-3.00079067694e-11_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,31,3.91175568865e-12_NTREAL)
    CALL SetChebyshevCoefficient(polynomial,32,-2.21155654398e-13_NTREAL)

    CALL FactorizedChebyshevCompute(ScaledMat, OutputMat, polynomial, &
         & f_sub_solver_parameters)

    !! Scale Back
    CALL ScaleDistributedSparseMatrix(OutputMat, &
         & REAL(2**(sigma_counter-1),NTREAL))

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructChebyshevPolynomial(polynomial)
    CALL DestructDistributedSparseMatrix(ScaledMat)
    CALL DestructDistributedSparseMatrix(IdentityMat)
    CALL DestructDistributedSparseMatrix(TempMat)
  END SUBROUTINE ComputeLogarithm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the logarithm of a matrix using a taylor series expansion.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = log(InputMat)
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeLogarithmTaylor(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: OutputMat
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: ScaledMat
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedSparseMatrix_t) :: Ak
    TYPE(DistributedSparseMatrix_t) :: IdentityMat
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Variables
    TYPE(IterativeSolverParameters_t) :: sub_solver_parameters
    REAL(NTREAL) :: e_min, e_max, spectral_radius
    REAL(NTREAL) :: sigma_val
    REAL(NTREAL) :: taylor_denom
    INTEGER :: sigma_counter
    INTEGER :: counter

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF
    CALL ConvertFixedToIterative(solver_parameters, sub_solver_parameters)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Logarithm Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Taylor")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL GershgorinBounds(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0
    sigma_counter = 1
    CALL CopyDistributedSparseMatrix(InputMat,ScaledMat)
    !do while (spectral_radius/sigma_val .gt. 1.1e-5)
    DO WHILE (spectral_radius/sigma_val .GT. 1.1e-7)
       CALL SquareRoot(ScaledMat,TempMat,sub_solver_parameters)
       CALL CopyDistributedSparseMatrix(TempMat,ScaledMat)
       CALL GershgorinBounds(ScaledMat, e_min, e_max)
       spectral_radius = MAX(ABS(e_min), ABS(e_max))
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL ConstructEmptyDistributedSparseMatrix(IdentityMat, &
         & InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(IdentityMat)

    !! Setup Matrices
    CALL IncrementDistributedSparseMatrix(IdentityMat,ScaledMat, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL CopyDistributedSparseMatrix(IdentityMat,Ak)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(ScaledMat, ScaledMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Ak, Ak, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Expand taylor series.
    CALL CopyDistributedSparseMatrix(ScaledMat,OutputMat)
    DO counter=2,10
       IF (MOD(counter,2) .EQ. 0) THEN
          taylor_denom = -1 * counter
       ELSE
          taylor_denom = counter
       END IF
       CALL DistributedGemm(Ak,ScaledMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparseMatrix(TempMat,Ak)
       CALL IncrementDistributedSparseMatrix(Ak, OutputMat, &
            & alpha_in=1.0/taylor_denom)
    END DO

    !! Undo scaling.
    CALL ScaleDistributedSparseMatrix(OutputMat,REAL(2**sigma_counter,NTREAL))

    !! Undo load balancing.
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructDistributedSparseMatrix(ScaledMat)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(IdentityMat)
    CALL DestructDistributedSparseMatrix(Ak)
  END SUBROUTINE ComputeLogarithmTaylor
END MODULE ExponentialSolversModule
