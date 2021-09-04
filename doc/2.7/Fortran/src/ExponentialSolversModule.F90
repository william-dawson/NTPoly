!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Matrix Exponentials and Logarithms.
MODULE ExponentialSolversModule
  USE ChebyshevSolversModule, ONLY : ChebyshevPolynomial_t, Compute, &
       & ConstructPolynomial, DestructPolynomial, FactorizedCompute, &
       & SetCoefficient
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds, PowerBounds
  USE LinearSolversModule, ONLY : CGSolver
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteElement
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, ScaleMatrix, &
       & IncrementMatrix
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, FillMatrixIdentity, PrintMatrixInformation
  USE RootSolversModule, ONLY : ComputeRoot
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  USE SquareRootSolversModule, ONLY : SquareRoot
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
  SUBROUTINE ComputeExponential(InputMat, OutputMat, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = exp(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    TYPE(SolverParameters_t) :: sub_solver_parameters
    TYPE(SolverParameters_t) :: psub_solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: ScaledMat
    TYPE(Matrix_ps) :: TempMat
    TYPE(MatrixMemoryPool_p) :: pool
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
       solver_parameters = SolverParameters_t()
    END IF
    sub_solver_parameters = solver_parameters
    psub_solver_parameters = solver_parameters
    psub_solver_parameters%max_iterations = 10

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Chebyshev")
       CALL PrintParameters(solver_parameters)
    END IF

    CALL ConstructEmptyMatrix(OutputMat, InputMat)

    !! Scale the matrix
    CALL PowerBounds(InputMat,spectral_radius,psub_solver_parameters)
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 1.0)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO
    CALL CopyMatrix(InputMat, ScaledMat)
    CALL ScaleMatrix(ScaledMat,1.0/sigma_val)
    sub_solver_parameters%threshold = sub_solver_parameters%threshold/sigma_val

    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Sigma", VALUE=sigma_val)
    END IF

    !! Expand Chebyshev Series
    CALL ConstructPolynomial(polynomial,16)
    CALL SetCoefficient(polynomial,1,1.266065877752007e+00_NTREAL)
    CALL SetCoefficient(polynomial,2,1.130318207984970e+00_NTREAL)
    CALL SetCoefficient(polynomial,3,2.714953395340771e-01_NTREAL)
    CALL SetCoefficient(polynomial,4,4.433684984866504e-02_NTREAL)
    CALL SetCoefficient(polynomial,5,5.474240442092110e-03_NTREAL)
    CALL SetCoefficient(polynomial,6,5.429263119148932e-04_NTREAL)
    CALL SetCoefficient(polynomial,7,4.497732295351912e-05_NTREAL)
    CALL SetCoefficient(polynomial,8,3.198436462630565e-06_NTREAL)
    CALL SetCoefficient(polynomial,9,1.992124801999838e-07_NTREAL)
    CALL SetCoefficient(polynomial,10,1.103677287249654e-08_NTREAL)
    CALL SetCoefficient(polynomial,11,5.505891628277851e-10_NTREAL)
    CALL SetCoefficient(polynomial,12,2.498021534339559e-11_NTREAL)
    CALL SetCoefficient(polynomial,13,1.038827668772902e-12_NTREAL)
    CALL SetCoefficient(polynomial,14,4.032447357431817e-14_NTREAL)
    CALL SetCoefficient(polynomial,15,2.127980007794583e-15_NTREAL)
    CALL SetCoefficient(polynomial,16,-1.629151584468762e-16_NTREAL)

    CALL Compute(ScaledMat,OutputMat,polynomial,sub_solver_parameters)
    !CALL FactorizedChebyshevCompute(ScaledMat,OutputMat,polynomial, &
    !     & sub_solver_parameters)

    !! Undo the scaling by squaring at the end.
    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    DO counter=1,sigma_counter-1
       CALL MatrixMultiply(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,OutputMat)
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
    CALL DestructPolynomial(polynomial)
    CALL DestructMatrix(ScaledMat)
    CALL DestructMatrix(TempMat)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE ComputeExponential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix using a pade approximation.
  !> Be warned, the pade method can result in a lot of intermediate fill.
  SUBROUTINE ComputeExponentialPade(InputMat, OutputMat, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = exp(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    TYPE(SolverParameters_t) :: sub_solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: ScaledMat
    TYPE(Matrix_ps) :: IdentityMat
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: B1, B2, B3
    TYPE(Matrix_ps) :: P1, P2
    TYPE(Matrix_ps) :: LeftMat, RightMat
    TYPE(MatrixMemoryPool_p) :: pool
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
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Pade")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Setup
    CALL ConstructEmptyMatrix(IdentityMat, InputMat)
    CALL FillMatrixIdentity(IdentityMat)

    !! Scale the matrix
    spectral_radius = MatrixNorm(InputMat)
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 1.0)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO
    CALL CopyMatrix(InputMat, ScaledMat)
    CALL ScaleMatrix(ScaledMat,1.0/sigma_val)
    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Sigma", VALUE=sigma_val)
       CALL WriteElement(key="Scaling_Steps", VALUE=sigma_counter)
    END IF

    !! Sub Solver Parameters
    sub_solver_parameters = solver_parameters
    sub_solver_parameters%threshold = sub_solver_parameters%threshold/sigma_val

    !! Power Matrices
    CALL MatrixMultiply(ScaledMat, ScaledMat, B1, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(B1, B1, B2, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(B2, B2, B3, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)

    !! Polynomials - 1
    CALL CopyMatrix(IdentityMat, P1)
    CALL ScaleMatrix(P1,17297280.0_NTREAL)
    CALL IncrementMatrix(B1, P1, alpha_in=1995840.0_NTREAL)
    CALL IncrementMatrix(B2, P1, alpha_in=25200.0_NTREAL)
    CALL IncrementMatrix(B3, P1, alpha_in=56.0_NTREAL)
    !! Polynomials - 2
    CALL CopyMatrix(IdentityMat, TempMat)
    CALL ScaleMatrix(TempMat,8648640.0_NTREAL)
    CALL IncrementMatrix(B1, TempMat, alpha_in=277200.0_NTREAL)
    CALL IncrementMatrix(B2, TempMat, alpha_in=1512.0_NTREAL)
    CALL IncrementMatrix(B3, TempMat)
    CALL MatrixMultiply(ScaledMat, TempMat, P2, &
         & threshold_in=sub_solver_parameters%threshold, memory_pool_in=pool)

    !! Left and Right
    CALL CopyMatrix(P1, LeftMat)
    CALL IncrementMatrix(P2, LeftMat, -1.0_NTREAL)
    CALL CopyMatrix(P1, RightMat)
    CALL IncrementMatrix(P2, RightMat, 1.0_NTREAL)

    CALL CGSolver(LeftMat, OutputMat, RightMat, sub_solver_parameters)

    !! Undo the scaling by squaring at the end.
    DO counter=1,sigma_counter-1
       CALL MatrixMultiply(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,OutputMat)
    END DO

    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(ScaledMat)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(B1)
    CALL DestructMatrix(B2)
    CALL DestructMatrix(B3)
    CALL DestructMatrix(P1)
    CALL DestructMatrix(P2)
    CALL DestructMatrix(LeftMat)
    CALL DestructMatrix(RightMat)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE ComputeExponentialPade
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix using a taylor series expansion.
  !> This is only really useful if you have a very small spectrum, because
  !> quite a bit of scaling is required.
  SUBROUTINE ComputeExponentialTaylor(InputMat, OutputMat, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = exp(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    TYPE(SolverParameters_t) :: psub_solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: ScaledMat
    TYPE(Matrix_ps) :: Ak
    TYPE(Matrix_ps) :: TempMat
    TYPE(MatrixMemoryPool_p) :: pool
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
       solver_parameters = SolverParameters_t()
    END IF
    psub_solver_parameters = solver_parameters
    psub_solver_parameters%max_iterations = 10

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Exponential Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Taylor")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL PowerBounds(InputMat,spectral_radius,psub_solver_parameters)

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 3.0e-8)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL CopyMatrix(InputMat, ScaledMat)
    CALL ScaleMatrix(ScaledMat,1.0/sigma_val)

    CALL ConstructEmptyMatrix(OutputMat, InputMat)
    CALL FillMatrixIdentity(OutputMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(ScaledMat, ScaledMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Expand Taylor Series
    taylor_denom = 1.0
    CALL CopyMatrix(OutputMat, Ak)
    DO counter=1,10
       taylor_denom = taylor_denom * counter
       CALL MatrixMultiply(Ak,ScaledMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,Ak)
       CALL IncrementMatrix(Ak,OutputMat)
    END DO

    DO counter=1,sigma_counter-1
       CALL MatrixMultiply(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,OutputMat)
    END DO

    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(ScaledMat)
    CALL DestructMatrix(Ak)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE ComputeExponentialTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the logarithm of a matrix.
  SUBROUTINE ComputeLogarithm(InputMat, OutputMat, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = exp(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: ScaledMat
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: IdentityMat
    !! For Chebyshev Expansion
    TYPE(ChebyshevPolynomial_t) :: polynomial
    !! Local Variables
    TYPE(SolverParameters_t) :: i_sub_solver_parameters
    TYPE(SolverParameters_t) :: p_sub_solver_parameters
    TYPE(SolverParameters_t) :: f_sub_solver_parameters
    REAL(NTREAL) :: spectral_radius
    INTEGER :: sigma_val
    INTEGER :: sigma_counter

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF
    i_sub_solver_parameters = solver_parameters
    p_sub_solver_parameters = solver_parameters
    p_sub_solver_parameters%max_iterations=16
    f_sub_solver_parameters = solver_parameters

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Logarithm Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Chebyshev")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Setup
    CALL ConstructEmptyMatrix(IdentityMat, InputMat)
    CALL FillMatrixIdentity(IdentityMat)

    !! Copy to a temporary matrix for scaling.
    CALL CopyMatrix(InputMat,ScaledMat)

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
       CALL WriteElement(key="Sigma", VALUE=sigma_val)
    END IF
    f_sub_solver_parameters%threshold = &
         & f_sub_solver_parameters%threshold/REAL(2**(sigma_counter-1),NTREAL)
    CALL ComputeRoot(InputMat, ScaledMat, sigma_val, i_sub_solver_parameters)

    !! Shift Scaled Matrix
    CALL IncrementMatrix(IdentityMat,ScaledMat, &
         & alpha_in=REAL(-1.0,NTREAL))

    !! Expand Chebyshev Series
    CALL ConstructPolynomial(polynomial,32)
    CALL SetCoefficient(polynomial,1,-0.485101351704_NTREAL)
    CALL SetCoefficient(polynomial,2,1.58828112379_NTREAL)
    CALL SetCoefficient(polynomial,3,-0.600947731795_NTREAL)
    CALL SetCoefficient(polynomial,4,0.287304748177_NTREAL)
    CALL SetCoefficient(polynomial,5,-0.145496447103_NTREAL)
    CALL SetCoefficient(polynomial,6,0.0734013668818_NTREAL)
    CALL SetCoefficient(polynomial,7,-0.0356277942958_NTREAL)
    CALL SetCoefficient(polynomial,8,0.0161605505166_NTREAL)
    CALL SetCoefficient(polynomial,9,-0.0066133591188_NTREAL)
    CALL SetCoefficient(polynomial,10,0.00229833505456_NTREAL)
    CALL SetCoefficient(polynomial,11,-0.000577804103964_NTREAL)
    CALL SetCoefficient(polynomial,12,2.2849332964e-05_NTREAL)
    CALL SetCoefficient(polynomial,13,8.37426826403e-05_NTREAL)
    CALL SetCoefficient(polynomial,14,-6.10822859027e-05_NTREAL)
    CALL SetCoefficient(polynomial,15,2.58132364523e-05_NTREAL)
    CALL SetCoefficient(polynomial,16,-5.87577322647e-06_NTREAL)
    CALL SetCoefficient(polynomial,17,-8.56711062722e-07_NTREAL)
    CALL SetCoefficient(polynomial,18,1.52066488969e-06_NTREAL)
    CALL SetCoefficient(polynomial,19,-7.12760496253e-07_NTREAL)
    CALL SetCoefficient(polynomial,20,1.23102245249e-07_NTREAL)
    CALL SetCoefficient(polynomial,21,6.03168259043e-08_NTREAL)
    CALL SetCoefficient(polynomial,22,-5.1865499826e-08_NTREAL)
    CALL SetCoefficient(polynomial,23,1.43185107512e-08_NTREAL)
    CALL SetCoefficient(polynomial,24,2.58449717089e-09_NTREAL)
    CALL SetCoefficient(polynomial,25,-3.73189861771e-09_NTREAL)
    CALL SetCoefficient(polynomial,26,1.18469334815e-09_NTREAL)
    CALL SetCoefficient(polynomial,27,1.51569931066e-10_NTREAL)
    CALL SetCoefficient(polynomial,28,-2.89595999673e-10_NTREAL)
    CALL SetCoefficient(polynomial,29,1.26720668874e-10_NTREAL)
    CALL SetCoefficient(polynomial,30,-3.00079067694e-11_NTREAL)
    CALL SetCoefficient(polynomial,31,3.91175568865e-12_NTREAL)
    CALL SetCoefficient(polynomial,32,-2.21155654398e-13_NTREAL)

    CALL FactorizedCompute(ScaledMat, OutputMat, polynomial, &
         & f_sub_solver_parameters)

    !! Scale Back
    CALL ScaleMatrix(OutputMat, &
         & REAL(2**(sigma_counter-1),NTREAL))

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructPolynomial(polynomial)
    CALL DestructMatrix(ScaledMat)
    CALL DestructMatrix(IdentityMat)
    CALL DestructMatrix(TempMat)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE ComputeLogarithm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the logarithm of a matrix using a taylor series expansion.
  SUBROUTINE ComputeLogarithmTaylor(InputMat, OutputMat, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = exp(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: ScaledMat
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: Ak
    TYPE(Matrix_ps) :: IdentityMat
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Variables
    TYPE(SolverParameters_t) :: sub_solver_parameters
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
       solver_parameters = SolverParameters_t()
    END IF
    sub_solver_parameters = solver_parameters

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Logarithm Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Taylor")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL GershgorinBounds(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0
    sigma_counter = 1
    CALL CopyMatrix(InputMat,ScaledMat)
    !do while (spectral_radius/sigma_val .gt. 1.1e-5)
    DO WHILE (spectral_radius/sigma_val .GT. 1.1e-7)
       CALL SquareRoot(ScaledMat,TempMat,sub_solver_parameters)
       CALL CopyMatrix(TempMat,ScaledMat)
       CALL GershgorinBounds(ScaledMat, e_min, e_max)
       spectral_radius = MAX(ABS(e_min), ABS(e_max))
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL ConstructEmptyMatrix(IdentityMat, InputMat)
    CALL FillMatrixIdentity(IdentityMat)

    !! Setup Matrices
    CALL IncrementMatrix(IdentityMat,ScaledMat, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL CopyMatrix(IdentityMat,Ak)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(ScaledMat, ScaledMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Ak, Ak, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Expand taylor series.
    CALL CopyMatrix(ScaledMat,OutputMat)
    DO counter=2,10
       IF (MOD(counter,2) .EQ. 0) THEN
          taylor_denom = -1 * counter
       ELSE
          taylor_denom = counter
       END IF
       CALL MatrixMultiply(Ak,ScaledMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,Ak)
       CALL IncrementMatrix(Ak, OutputMat, &
            & alpha_in=1.0/taylor_denom)
    END DO

    !! Undo scaling.
    CALL ScaleMatrix(OutputMat,REAL(2**sigma_counter,NTREAL))

    !! Undo load balancing.
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(ScaledMat)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(IdentityMat)
    CALL DestructMatrix(Ak)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE ComputeLogarithmTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ExponentialSolversModule
