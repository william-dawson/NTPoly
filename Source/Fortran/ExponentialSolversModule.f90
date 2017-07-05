!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Matrix Exponentials and Logarithms.
MODULE ExponentialSolversModule
  USE ChebyshevSolversModule
  USE DataTypesModule
#ifdef NOBLOCK
  USE DistributedSparseMatrixModule
#else
  USE DistributedBlockedSparseMatrixModule
#endif
  USE DistributedMatrixMemoryPoolModule
  USE FixedSolversModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ParameterConverterModule
  USE ProcessGridModule
  USE RootSolversModule
  USE SquareRootSolversModule
  USE TimerModule
  USE mpi
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

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Chebyshev")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    CALL ConstructEmpty(OutputMat, InputMat%actual_matrix_dimension)

    !! Scale the matrix
    CALL EigenCircle(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))
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
    END IF

    !! Expand Chebyshev Series
    CALL ConstructChebyshevPolynomial(polynomial,16)
    CALL SetChebyshevCoefficient(polynomial,1, &
         & REAL(1.266065877752007e+00,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,2, &
         & REAL(1.130318207984970e+00,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,3, &
         & REAL(2.714953395340764e-01,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,4, &
         & REAL(4.433684984866457e-02,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,5, &
         & REAL(5.474240442093185e-03,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,6, &
         & REAL(5.429263119148893e-04,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,7, &
         & REAL(4.497732295384419e-05,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,8, &
         & REAL(3.198436462646902e-06,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,9, &
         & REAL(1.992124791702536e-07,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,10, &
         & REAL(1.103677179722673e-08,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,11, &
         & REAL(5.505894210845503e-10,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,12, &
         & REAL(2.497953174170452e-11,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,13, &
         & REAL(1.038796461263896e-12,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,14, &
         & REAL(3.933483178134044e-14,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,15, &
         & REAL(4.194318048730632e-16,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,16, &
         & REAL(-4.462345204611966e-16,NTREAL))

    !CALL ChebyshevCompute(ScaledMat,OutputMat,polynomial,solver_parameters)
    CALL FactorizedChebyshevCompute(ScaledMat,OutputMat,polynomial, &
         & solver_parameters)

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
    CALL EigenCircle(InputMat, e_min, e_max)
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
    CALL EigenCircle(ScaledMat, e_min, e_max)
    sigma_val = 1
    sigma_counter = 1
    spectral_radius = MAX(ABS(e_min), ABS(e_max))
    DO WHILE (spectral_radius .GT. 4)
       spectral_radius = SQRT(spectral_radius)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO
    CALL ComputeRoot(InputMat, ScaledMat, sigma_val, i_sub_solver_parameters)

    !! Shift Scaled Matrix
    CALL IncrementDistributedSparseMatrix(IdentityMat,ScaledMat, &
         & alpha_in=REAL(-3.0,NTREAL))

    !! Expand Chebyshev Series
    CALL ConstructChebyshevPolynomial(polynomial,16)
    CALL SetChebyshevCoefficient(polynomial,1, &
         & REAL(1.069599993479196e+00,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,2, &
         & REAL(3.431457505076075e-01,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,3, &
         & REAL(-2.943725152284926e-02,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,4, &
         & REAL(3.367089255551476e-03,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,5, &
         & REAL(-4.332758885990646e-04,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,6, &
         & REAL(5.947071197674775e-05,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,7, &
         & REAL(-8.502967529983385e-06,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,8, &
         & REAL(1.250467349906989e-06,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,9, &
         & REAL(-1.877279855374288e-07,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,10, &
         & REAL(2.863023858749194e-08,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,11, &
         & REAL(-4.420947514783171e-09,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,12, &
         & REAL(6.895489131189299e-10,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,13, &
         & REAL(-1.084419051323420e-10,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,14, &
         & REAL(1.716514621065752e-11,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,15, &
         & REAL(-2.729866800230529e-12,NTREAL))
    CALL SetChebyshevCoefficient(polynomial,16, &
         & REAL(4.277875177870295e-13,NTREAL))

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
    CALL EigenCircle(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0
    sigma_counter = 1
    CALL CopyDistributedSparseMatrix(InputMat,ScaledMat)
    !do while (spectral_radius/sigma_val .gt. 1.1e-5)
    DO WHILE (spectral_radius/sigma_val .GT. 1.1e-7)
       CALL SquareRoot(ScaledMat,TempMat,sub_solver_parameters)
       CALL CopyDistributedSparseMatrix(TempMat,ScaledMat)
       CALL EigenCircle(ScaledMat, e_min, e_max)
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
