!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Matrix Exponentials and Logarithms.
MODULE ExponentialSolversModule
  USE ChebyshevSolversModule
  USE DataTypesModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixModule
  USE EigenBoundsModule
  USE FixedSolversModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ParameterConverterModule
  USE ProcessGridModule
  USE RootSolversModule
  USE SquareRootSolversModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: ComputeExponential
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
    TYPE(DistributedSparseMatrix), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(inout) :: OutputMat
    TYPE(FixedSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters) :: solver_parameters
    TYPE(FixedSolverParameters) :: sub_solver_parameters
    TYPE(IterativeSolverParameters) :: i_sub_solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: ScaledMat
    TYPE(DistributedSparseMatrix) :: TempMat
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! For Chebyshev Expansion
    TYPE(ChebyshevPolynomial_t) :: polynomial
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max, spectral_radius
    REAL(NTREAL) :: sigma_val
    INTEGER :: sigma_counter
    INTEGER :: counter
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters()
    END IF
    sub_solver_parameters = solver_parameters
    CALL ConvertFixedToIterative(solver_parameters, i_sub_solver_parameters)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Chebyshev")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    CALL ConstructEmpty(OutputMat, InputMat%actual_matrix_dimension)

    !! Scale the matrix
    ! CALL GershgorinBounds(InputMat, e_min, e_max)
    ! spectral_radius = MAX(ABS(e_min), ABS(e_max))
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
    CALL SetChebyshevCoefficient(polynomial,1, &
         & REAL(1.266065877752007e+00_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,2, &
         & REAL(1.130318207984970e+00_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,3, &
         & REAL(2.714953395340771e-01_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,4, &
         & REAL(4.433684984866504e-02_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,5, &
         & REAL(5.474240442092110e-03_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,6, &
         & REAL(5.429263119148932e-04_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,7, &
         & REAL(4.497732295351912e-05_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,8, &
         & REAL(3.198436462630565e-06_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,9, &
         & REAL(1.992124801999838e-07_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,10, &
         & REAL(1.103677287249654e-08_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,11, &
         & REAL(5.505891628277851e-10_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,12, &
         & REAL(2.498021534339559e-11_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,13, &
         & REAL(1.038827668772902e-12_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,14, &
         & REAL(4.032447357431817e-14_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,15, &
         & REAL(2.127980007794583e-15_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,16, &
         & REAL(-1.629151584468762e-16_16,NTREAL))

    !CALL ChebyshevCompute(ScaledMat,OutputMat,polynomial,solver_parameters)
    CALL FactorizedChebyshevCompute(ScaledMat,OutputMat,polynomial, &
         & sub_solver_parameters)

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
       CALL GetLoadBalance(OutputMat,min_size,max_size)
       sparsity = REAL(GetSize(OutputMat),KIND=NTREAL)/ &
            & (OutputMat%actual_matrix_dimension**2)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
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
  !> Compute the exponential of a matrix using a taylor series expansion.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = exp(InputMat)
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeExponentialTaylor(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(inout) :: OutputMat
    TYPE(FixedSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: ScaledMat
    TYPE(DistributedSparseMatrix) :: Ak
    TYPE(DistributedSparseMatrix) :: TempMat
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Variables
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
       solver_parameters = FixedSolverParameters()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
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
    DO WHILE (spectral_radius/sigma_val .GT. 3.0e-8)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL CopyDistributedSparseMatrix(InputMat, ScaledMat)
    CALL ScaleDistributedSparseMatrix(ScaledMat,1.0/sigma_val)

    CALL ConstructEmpty(OutputMat, InputMat%actual_matrix_dimension)
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
  !! @param[out] OutputMat = exp(InputMat)
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeLogarithm(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(IN)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(INOUT) :: OutputMat
    TYPE(FixedSolverParameters), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: ScaledMat
    TYPE(DistributedSparseMatrix) :: TempMat
    TYPE(DistributedSparseMatrix) :: IdentityMat
    !! For Chebyshev Expansion
    TYPE(ChebyshevPolynomial_t) :: polynomial
    !! Local Variables
    TYPE(IterativeSolverParameters) :: i_sub_solver_parameters
    TYPE(FixedSolverParameters) :: f_sub_solver_parameters
    REAL(NTREAL) :: e_min, e_max, spectral_radius
    INTEGER :: sigma_val
    INTEGER :: sigma_counter
    INTEGER :: target_root

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters()
    END IF
    CALL ConvertFixedToIterative(solver_parameters, i_sub_solver_parameters)
    f_sub_solver_parameters = solver_parameters

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Logarithm Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Chebyshev")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Setup
    CALL ConstructEmpty(IdentityMat, InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(IdentityMat)

    !! Copy to a temporary matrix for scaling.
    CALL CopyDistributedSparseMatrix(InputMat,ScaledMat)

    !! Compute The Scaling Factor
    CALL GershgorinBounds(ScaledMat, e_min, e_max)
    sigma_val = 1
    sigma_counter = 1
    spectral_radius = MAX(ABS(e_min), ABS(e_max))
    DO WHILE (spectral_radius .GT. SQRT(2.0))
       spectral_radius = SQRT(spectral_radius)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Sigma", int_value_in=sigma_val)
    END IF
    CALL ComputeRoot(InputMat, ScaledMat, sigma_val, i_sub_solver_parameters)

    !! Shift Scaled Matrix
    CALL IncrementDistributedSparseMatrix(IdentityMat,ScaledMat, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL GershgorinBounds(ScaledMat, e_min, e_max)

    !! Expand Chebyshev Series
    CALL ConstructChebyshevPolynomial(polynomial,32)
    CALL SetChebyshevCoefficient(polynomial,1, &
         & REAL(-0.485101351704_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,2, &
         & REAL(1.58828112379_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,3, &
         & REAL(-0.600947731795_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,4, &
         & REAL(0.287304748177_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,5, &
         & REAL(-0.145496447103_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,6, &
         & REAL(0.0734013668818_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,7, &
         & REAL(-0.0356277942958_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,8, &
         & REAL(0.0161605505166_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,9, &
         & REAL(-0.0066133591188_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,10, &
         & REAL(0.00229833505456_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,11, &
         & REAL(-0.000577804103964_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,12, &
         & REAL(2.2849332964e-05_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,13, &
         & REAL(8.37426826403e-05_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,14, &
         & REAL(-6.10822859027e-05_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,15, &
         & REAL(2.58132364523e-05_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,16, &
         & REAL(-5.87577322647e-06_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,17, &
         & REAL(-8.56711062722e-07_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,18, &
         & REAL(1.52066488969e-06_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,19, &
         & REAL(-7.12760496253e-07_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,20, &
         & REAL(1.23102245249e-07_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,21, &
         & REAL(6.03168259043e-08_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,22, &
         & REAL(-5.1865499826e-08_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,23, &
         & REAL(1.43185107512e-08_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,24, &
         & REAL(2.58449717089e-09_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,25, &
         & REAL(-3.73189861771e-09_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,26, &
         & REAL(1.18469334815e-09_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,27, &
         & REAL(1.51569931066e-10_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,28, &
         & REAL(-2.89595999673e-10_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,29, &
         & REAL(1.26720668874e-10_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,30, &
         & REAL(-3.00079067694e-11_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,31, &
         & REAL(3.91175568865e-12_16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,32, &
         & REAL(-2.21155654398e-13_16,NTREAL))

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
    TYPE(DistributedSparseMatrix), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(inout) :: OutputMat
    TYPE(FixedSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: ScaledMat
    TYPE(DistributedSparseMatrix) :: TempMat
    TYPE(DistributedSparseMatrix) :: Ak
    TYPE(DistributedSparseMatrix) :: IdentityMat
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Variables
    TYPE(IterativeSolverParameters) :: sub_solver_parameters
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
       solver_parameters = FixedSolverParameters()
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

    CALL ConstructEmpty(IdentityMat, InputMat%actual_matrix_dimension)
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
