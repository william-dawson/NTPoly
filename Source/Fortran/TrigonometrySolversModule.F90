!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Trigonometric functions of a Matrix.
MODULE TrigonometrySolversModule
  USE DataTypesModule
  USE MatrixMemoryPoolPModule
  USE MatrixPSAlgebraModule
  USE MatrixPSModule
  USE EigenBoundsModule
  USE FixedSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: Sine
  PUBLIC :: Cosine
  PUBLIC :: ScaleSquareTrigonometryTaylor
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the sine of a matrix.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] solver_parameters_in parameters for the solver, optional.
  SUBROUTINE Sine(InputMat, OutputMat, solver_parameters_in)
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(FixedSolverParameters_t),INTENT(in),OPTIONAL :: solver_parameters_in
    !! A temporary matrix to hold the transformation from sine to cosine.
    TYPE(Matrix_ps) :: ShiftedMat
    TYPE(Matrix_ps) :: IdentityMat
    REAL(NTREAL), PARAMETER :: PI = 4*ATAN(1.0)

    !! Shift
    CALL CopyMatrix(InputMat,ShiftedMat)
    CALL ConstructEmptyMatrix(IdentityMat, &
         InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL FillMatrixIdentity(IdentityMat)
    CALL IncrementMatrix(IdentityMat,ShiftedMat, &
         & alpha_in=REAL(-1.0*PI/2.0,NTREAL))
    CALL DestructMatrix(IdentityMat)

    IF (PRESENT(solver_parameters_in)) THEN
       CALL ScaleSquareTrigonometry(ShiftedMat, OutputMat, solver_parameters_in)
    ELSE
       CALL ScaleSquareTrigonometry(ShiftedMat, OutputMat)
    END IF
    CALL DestructMatrix(ShiftedMat)
  END SUBROUTINE Sine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the cosine of a matrix.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] solver_parameters_in parameters for the solver, optional.
  SUBROUTINE Cosine(InputMat, OutputMat, solver_parameters_in)
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(FixedSolverParameters_t),INTENT(in),OPTIONAL :: solver_parameters_in
    IF (PRESENT(solver_parameters_in)) THEN
       CALL ScaleSquareTrigonometry(InputMat, OutputMat, solver_parameters_in)
    ELSE
       CALL ScaleSquareTrigonometry(InputMat, OutputMat)
    END IF
  END SUBROUTINE Cosine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute trigonometric functions of a matrix using a taylor series.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE ScaleSquareTrigonometryTaylor(InputMat, OutputMat, &
       & solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(FixedSolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
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
       solver_parameters = FixedSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Trigonometry Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Taylor")
       CALL WriteCitation("higham2003computing")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL GershgorinBounds(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 3.0e-3)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL CopyMatrix(InputMat, ScaledMat)
    CALL ScaleMatrix(ScaledMat,1.0/sigma_val)
    CALL ConstructEmptyMatrix(OutputMat, &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL FillMatrixIdentity(OutputMat)
    CALL ConstructEmptyMatrix(IdentityMat, &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
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
    taylor_denom = -2.0
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
            & alpha_in=REAL(1.0/taylor_denom,NTREAL))
       taylor_denom = taylor_denom * (counter+1)
       taylor_denom = -1.0*taylor_denom*(counter+1)
    END DO

    !! Undo scaling
    DO counter=1,sigma_counter-1
       CALL MatrixMultiply(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,OutputMat)
       CALL ScaleMatrix(OutputMat,REAL(2.0,NTREAL))
       CALL IncrementMatrix(IdentityMat,OutputMat, &
            & REAL(-1.0,NTREAL))
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
  END SUBROUTINE ScaleSquareTrigonometryTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute trigonometric functions of a matrix.
  !! This method uses Chebyshev polynomials.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] compute_cosine_in cosine or sign (optional, default true).
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE ScaleSquareTrigonometry(InputMat, OutputMat, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(FixedSolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
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
       solver_parameters = FixedSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Trigonometry Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Chebyshev")
       CALL WriteCitation("serbin1980algorithm higham2003computing yau1993reducing")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL GershgorinBounds(InputMat, e_min, e_max)
    spectral_radius = MAX(ABS(e_min), ABS(e_max))

    !! Figure out how much to scale the matrix.
    sigma_val = 1.0
    sigma_counter = 1
    DO WHILE (spectral_radius/sigma_val .GT. 1.0)
       sigma_val = sigma_val * 2
       sigma_counter = sigma_counter + 1
    END DO

    CALL CopyMatrix(InputMat, ScaledMat)
    CALL ScaleMatrix(ScaledMat,1.0/sigma_val)
    CALL ConstructEmptyMatrix(OutputMat, &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL ConstructEmptyMatrix(IdentityMat, &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL FillMatrixIdentity(IdentityMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(ScaledMat, ScaledMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IdentityMat, IdentityMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Expand the Chebyshev Polynomial.
    coefficients(1) = 7.651976865579664e-01_16
    coefficients(2) = 0
    coefficients(3) = -2.298069698638004e-01_16
    coefficients(4) = 0
    coefficients(5) = 4.953277928219409e-03_16
    coefficients(6) = 0
    coefficients(7) = -4.187667600472235e-05_16
    coefficients(8) = 0
    coefficients(9) = 1.884468822397086e-07_16
    coefficients(10) = 0
    coefficients(11) = -5.261224549346905e-10_16
    coefficients(12) = 0
    coefficients(13) = 9.999906645345580e-13_16
    coefficients(14) = 0
    coefficients(15) = -2.083597362700025e-15_16
    coefficients(16) = 0
    coefficients(17) = 9.181480886537484e-17_16

    !! Basic T Values.
    CALL MatrixMultiply(ScaledMat,ScaledMat,T2,alpha_in=REAL(2.0,NTREAL),&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(IdentityMat,T2, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL MatrixMultiply(T2,T2,T4,alpha_in=REAL(2.0,NTREAL),&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(IdentityMat,T4, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL MatrixMultiply(T4,T2,T6,alpha_in=REAL(2.0,NTREAL),&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(T2,T6, &
         & alpha_in=REAL(-1.0,NTREAL))
    CALL MatrixMultiply(T6,T2,T8,alpha_in=REAL(2.0,NTREAL),&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL IncrementMatrix(T4,T8, &
         & alpha_in=REAL(-1.0,NTREAL))

    !! Contribution from the second half.
    CALL CopyMatrix(T8,OutputMat)
    CALL ScaleMatrix(OutputMat,0.5*coefficients(17))
    CALL IncrementMatrix(T6,OutputMat,&
         & alpha_in=REAL(0.5*coefficients(15),NTREAL))
    CALL IncrementMatrix(T4,OutputMat,&
         & alpha_in=REAL(0.5*coefficients(13),NTREAL))
    CALL IncrementMatrix(T2,OutputMat,&
         & alpha_in=REAL(0.5*coefficients(11),NTREAL))
    CALL MatrixMultiply(T8,OutputMat,TempMat,&
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Contribution from the first half.
    CALL CopyMatrix(T8,OutputMat)
    CALL ScaleMatrix(OutputMat,coefficients(9))
    CALL IncrementMatrix(T6,OutputMat,&
         & alpha_in=REAL(coefficients(7)+0.5*coefficients(11),NTREAL))
    CALL IncrementMatrix(T4,OutputMat,&
         & alpha_in=REAL(coefficients(5)+0.5*coefficients(13),NTREAL))
    CALL IncrementMatrix(T2,OutputMat,&
         & alpha_in=REAL(coefficients(3)+0.5*coefficients(15),NTREAL))
    CALL IncrementMatrix(IdentityMat,OutputMat,&
         & alpha_in=REAL(coefficients(1)+0.5*coefficients(17),NTREAL))

    CALL IncrementMatrix(TempMat,OutputMat)

    !! Undo scaling
    DO counter=1,sigma_counter-1
       CALL MatrixMultiply(OutputMat,OutputMat,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyMatrix(TempMat,OutputMat)
       CALL ScaleMatrix(OutputMat,REAL(2.0,NTREAL))
       CALL IncrementMatrix(IdentityMat,OutputMat, &
            & REAL(-1.0,NTREAL))
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
  END SUBROUTINE ScaleSquareTrigonometry
END MODULE TrigonometrySolversModule
