!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Trigonometric functions of a Matrix.
MODULE TrigonometrySolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteListElement, WriteElement
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, &
       & ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, FillMatrixIdentity
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: Sine
  PUBLIC :: Cosine
  PUBLIC :: ScaleSquareTrigonometryTaylor
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the sine of a matrix.
  SUBROUTINE Sine(InputMat, OutputMat, solver_parameters_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The resulting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t),INTENT(IN),OPTIONAL :: solver_parameters_in
    !! A temporary matrix to hold the transformation from sine to cosine.
    TYPE(Matrix_ps) :: ShiftedMat
    TYPE(Matrix_ps) :: IdentityMat
    REAL(NTREAL), PARAMETER :: PI = 4*ATAN(1.00_NTREAL)

    !! Shift
    CALL CopyMatrix(InputMat,ShiftedMat)
    CALL ConstructEmptyMatrix(IdentityMat, InputMat)
    CALL FillMatrixIdentity(IdentityMat)
    CALL IncrementMatrix(IdentityMat,ShiftedMat, &
         & alpha_in=REAL(-1.0_NTREAL*PI/2.0_NTREAL,NTREAL))
    CALL DestructMatrix(IdentityMat)

    IF (PRESENT(solver_parameters_in)) THEN
       CALL ScaleSquareTrigonometry(ShiftedMat, OutputMat, solver_parameters_in)
    ELSE
       CALL ScaleSquareTrigonometry(ShiftedMat, OutputMat)
    END IF

    !! Cleanup
    CALL DestructMatrix(ShiftedMat)
  END SUBROUTINE Sine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the cosine of a matrix.
  SUBROUTINE Cosine(InputMat, OutputMat, solver_parameters_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The resulting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t),INTENT(IN),OPTIONAL :: solver_parameters_in

    IF (PRESENT(solver_parameters_in)) THEN
       CALL ScaleSquareTrigonometry(InputMat, OutputMat, solver_parameters_in)
    ELSE
       CALL ScaleSquareTrigonometry(InputMat, OutputMat)
    END IF
  END SUBROUTINE Cosine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute trigonometric functions of a matrix using a taylor series.
  SUBROUTINE ScaleSquareTrigonometryTaylor(InputMat, OutputMat, &
       & solver_parameters_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The resulting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t),INTENT(IN),OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: ScaledMat
    TYPE(Matrix_ps) :: Ak
    TYPE(Matrix_ps) :: TempMat
    TYPE(MatrixMemoryPool_p) :: pool
    TYPE(Matrix_ps) :: IdentityMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max, spectral_radius
    REAL(NTREAL) :: sigma_val
    REAL(NTREAL) :: taylor_denom
    INTEGER :: sigma_counter
    INTEGER :: counter

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Trigonometry Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Taylor")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("higham2003computing")
       CALL ExitSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL GershgorinBounds(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0_NTREAL
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 3.0e-3_NTREAL)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL CopyMatrix(InputMat, ScaledMat)
    CALL ScaleMatrix(ScaledMat,1.0_NTREAL/sigma_val)
    CALL ConstructEmptyMatrix(OutputMat, InputMat)
    CALL FillMatrixIdentity(OutputMat)
    CALL ConstructEmptyMatrix(IdentityMat, InputMat)
    CALL FillMatrixIdentity(IdentityMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(ScaledMat, ScaledMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IdentityMat, IdentityMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Square the scaled matrix.
    taylor_denom = -2.0_NTREAL
    CALL CopyMatrix(OutputMat, Ak)
    CALL MatrixMultiply(ScaledMat,ScaledMat,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL CopyMatrix(TempMat,ScaledMat)

    !! Expand Taylor Series
    DO counter=2,40,2
       CALL MatrixMultiply(Ak,ScaledMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,Ak)
       CALL IncrementMatrix(Ak,OutputMat, &
            & alpha_in=REAL(1.0_NTREAL/taylor_denom,NTREAL))
       taylor_denom = taylor_denom * (counter+1)
       taylor_denom = -1.0_NTREAL*taylor_denom*(counter+1)
    END DO

    !! Undo scaling
    DO counter=1,sigma_counter-1
       CALL MatrixMultiply(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,OutputMat)
       CALL ScaleMatrix(OutputMat,REAL(2.0_NTREAL,NTREAL))
       CALL IncrementMatrix(IdentityMat,OutputMat, &
            & REAL(-1.0_NTREAL,NTREAL))
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
    CALL DestructMatrix(IdentityMat)
    CALL DestructMatrix(Ak)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE ScaleSquareTrigonometryTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute trigonometric functions of a matrix.
  !> This method uses Chebyshev polynomials.
  SUBROUTINE ScaleSquareTrigonometry(InputMat, OutputMat, solver_parameters_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The resulting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver.
    TYPE(SolverParameters_t),INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: ScaledMat
    TYPE(Matrix_ps) :: TempMat
    TYPE(MatrixMemoryPool_p) :: pool
    TYPE(Matrix_ps) :: IdentityMat
    !! For Chebyshev Expansion
    REAL(NTREAL), DIMENSION(17) :: coefficients
    TYPE(Matrix_ps) :: T2
    TYPE(Matrix_ps) :: T4
    TYPE(Matrix_ps) :: T6
    TYPE(Matrix_ps) :: T8
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max, spectral_radius
    REAL(NTREAL) :: sigma_val
    INTEGER :: sigma_counter
    INTEGER :: counter

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Trigonometry Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Chebyshev")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("serbin1980algorithm")
       CALL WriteListElement("higham2003computing")
       CALL WriteListElement("yau1993reducing")
       CALL ExitSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL GershgorinBounds(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0_NTREAL
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 1.0_NTREAL)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL CopyMatrix(InputMat, ScaledMat)
    CALL ScaleMatrix(ScaledMat,1.0_NTREAL/sigma_val)
    CALL ConstructEmptyMatrix(OutputMat, InputMat)
    CALL ConstructEmptyMatrix(IdentityMat, InputMat)
    CALL FillMatrixIdentity(IdentityMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(ScaledMat, ScaledMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IdentityMat, IdentityMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Expand the Chebyshev Polynomial.
    coefficients(1) = 7.651976865579664e-01_NTREAL
    coefficients(2) = 0_NTREAL
    coefficients(3) = -2.298069698638004e-01_NTREAL
    coefficients(4) = 0_NTREAL
    coefficients(5) = 4.953277928219409e-03_NTREAL
    coefficients(6) = 0_NTREAL
    coefficients(7) = -4.187667600472235e-05_NTREAL
    coefficients(8) = 0_NTREAL
    coefficients(9) = 1.884468822397086e-07_NTREAL
    coefficients(10) = 0_NTREAL
    coefficients(11) = -5.261224549346905e-10_NTREAL
    coefficients(12) = 0_NTREAL
    coefficients(13) = 9.999906645345580e-13_NTREAL
    coefficients(14) = 0_NTREAL
    coefficients(15) = -2.083597362700025e-15_NTREAL
    coefficients(16) = 0_NTREAL
    coefficients(17) = 9.181480886537484e-17_NTREAL

    !! Basic T Values.
    CALL MatrixMultiply(ScaledMat,ScaledMat,T2,alpha_in=2.0_NTREAL,&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(IdentityMat,T2, alpha_in=-1.0_NTREAL)
    CALL MatrixMultiply(T2,T2,T4,alpha_in=2.0_NTREAL,&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(IdentityMat,T4, alpha_in=-1.0_NTREAL)
    CALL MatrixMultiply(T4,T2,T6,alpha_in=2.0_NTREAL,&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(T2,T6, alpha_in=-1.0_NTREAL)
    CALL MatrixMultiply(T6,T2,T8,alpha_in=2.0_NTREAL,&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(T4,T8, alpha_in=-1.0_NTREAL)

    !! Contribution from the second half.
    CALL CopyMatrix(T8,OutputMat)
    CALL ScaleMatrix(OutputMat,0.5_NTREAL*coefficients(17))
    CALL IncrementMatrix(T6,OutputMat,alpha_in=0.5_NTREAL*coefficients(15))
    CALL IncrementMatrix(T4,OutputMat,alpha_in=0.5_NTREAL*coefficients(13))
    CALL IncrementMatrix(T2,OutputMat,alpha_in=0.5_NTREAL*coefficients(11))
    CALL MatrixMultiply(T8,OutputMat,TempMat,&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Contribution from the first half.
    CALL CopyMatrix(T8,OutputMat)
    CALL ScaleMatrix(OutputMat,coefficients(9))
    CALL IncrementMatrix(T6,OutputMat,&
         & alpha_in=coefficients(7)+0.5_NTREAL*coefficients(11))
    CALL IncrementMatrix(T4,OutputMat,&
         & alpha_in=coefficients(5)+0.5_NTREAL*coefficients(13))
    CALL IncrementMatrix(T2,OutputMat,&
         & alpha_in=coefficients(3)+0.5_NTREAL*coefficients(15))
    CALL IncrementMatrix(IdentityMat,OutputMat,&
         & alpha_in=coefficients(1)+0.5_NTREAL*coefficients(17))

    CALL IncrementMatrix(TempMat,OutputMat)

    !! Undo scaling
    DO counter=1,sigma_counter-1
       CALL MatrixMultiply(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,OutputMat)
       CALL ScaleMatrix(OutputMat,2.0_NTREAL)
       CALL IncrementMatrix(IdentityMat,OutputMat,-1.0_NTREAL)
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
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(IdentityMat)
    CALL DestructMatrix(T2)
    CALL DestructMatrix(T4)
    CALL DestructMatrix(T6)
    CALL DestructMatrix(T8)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE ScaleSquareTrigonometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TrigonometrySolversModule
