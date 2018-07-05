!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Matrix functions based on Chebyshev polynomials.
MODULE ChebyshevSolversModule
  USE DataTypesModule
  USE LoadBalancerModule
  USE LoggingModule
  USE PMatrixMemoryPoolModule
  USE PSMatrixAlgebraModule
  USE PSMatrixModule
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype that represents a Chebyshev polynomial.
  TYPE, PUBLIC :: ChebyshevPolynomial_t
     !> Coefficients of the polynomial.
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: coefficients
  END TYPE ChebyshevPolynomial_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Polynomial type
  PUBLIC :: ConstructChebyshevPolynomial
  PUBLIC :: DestructChebyshevPolynomial
  PUBLIC :: SetChebyshevCoefficient
  !! Solvers
  PUBLIC :: ChebyshevCompute
  PUBLIC :: FactorizedChebyshevCompute
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a Chebyshev polynomial object.
  !! @param[inout] this the polynomial to construct.
  !! @param[in] degree of the polynomial.
  PURE SUBROUTINE ConstructChebyshevPolynomial(this, degree)
    !! Parameters
    TYPE(ChebyshevPolynomial_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0
  END SUBROUTINE ConstructChebyshevPolynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  !! @param[inout] this the polynomial to destruct.
  PURE SUBROUTINE DestructChebyshevPolynomial(this)
    !! Parameters
    TYPE(ChebyshevPolynomial_t), INTENT(INOUT) :: this
    IF (ALLOCATED(this%coefficients)) THEN
       DEALLOCATE(this%coefficients)
    END IF
  END SUBROUTINE DestructChebyshevPolynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set a coefficient of a Chebyshev polynomial.
  !! @param[inout] this the polynomial to set.
  !! @param[in] degree for which to set the coefficient.
  !! @param[in] coefficient value.
  SUBROUTINE SetChebyshevCoefficient(this, degree, coefficient)
    !! Parameters
    TYPE(ChebyshevPolynomial_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: degree
    REAL(NTREAL), INTENT(IN) :: coefficient

    this%coefficients(degree) = coefficient
  END SUBROUTINE SetChebyshevCoefficient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Chebyshev Polynomial of the matrix.
  !! This method uses the standard Chebyshev Polynomial expansion.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in] poly the Chebyshev polynomial to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ChebyshevCompute(InputMat, OutputMat, poly, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    TYPE(ChebyshevPolynomial_t), INTENT(IN) :: poly
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(Matrix_ps) :: Tk
    TYPE(Matrix_ps) :: Tkminus1
    TYPE(Matrix_ps) :: Tkminus2
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Variables
    INTEGER :: degree
    INTEGER :: counter

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Chebyshev Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Standard")
       CALL WriteElement(key="Degree", int_value_in=degree)
       CALL PrintParameters(solver_parameters)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL CopyMatrix(InputMat,BalancedInput)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! First Term
    CALL CopyMatrix(Identity,Tkminus2)
    IF (degree == 1) THEN
       CALL CopyMatrix(Tkminus2,OutputMat)
       CALL ScaleMatrix(OutputMat,poly%coefficients(1))
    ELSE
       CALL CopyMatrix(BalancedInput,Tkminus1)
       CALL CopyMatrix(Tkminus2,OutputMat)
       CALL ScaleMatrix(OutputMat,poly%coefficients(1))
       CALL IncrementMatrix(Tkminus1,OutputMat, &
            & alpha_in=poly%coefficients(2))
       IF (degree > 2) THEN
          CALL MatrixMultiply(BalancedInput, Tkminus1, Tk, &
               & alpha_in=REAL(2.0,NTREAL), &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
          CALL IncrementMatrix(Tkminus2,Tk,REAL(-1.0,NTREAL))
          CALL IncrementMatrix(Tk, OutputMat, &
               & alpha_in=poly%coefficients(3))
          DO counter = 4, degree
             CALL CopyMatrix(Tkminus1,Tkminus2)
             CALL CopyMatrix(Tk,Tkminus1)
             CALL MatrixMultiply(BalancedInput, Tkminus1, Tk, &
                  & alpha_in=REAL(2.0,NTREAL), &
                  & threshold_in=solver_parameters%threshold, &
                  & memory_pool_in=pool)
             CALL IncrementMatrix(Tkminus2,Tk, &
                  & REAL(-1.0,NTREAL))
             CALL IncrementMatrix(Tk, OutputMat, &
                  & alpha_in=poly%coefficients(counter))
          END DO
       END IF
    END IF
    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(Tk)
    CALL DestructMatrix(Tkminus1)
    CALL DestructMatrix(Tkminus2)
    CALL DestructMatrix(BalancedInput)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE ChebyshevCompute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Chebyshev Polynomial of the matrix.
  !! This version first factors the Chebyshev Polynomial and computes the
  !! function using a divide and conquer algorithm. Based on a simplified
  !! version of the first method in \cite liang2003improved .
  !! @param[in]  InputMat the input matrix
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in]  poly the Chebyshev polynomial to compute.
  !! @param[in]  solver_parameters_in parameters for the solver (optional).
  SUBROUTINE FactorizedChebyshevCompute(InputMat, OutputMat, poly, &
       & solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    TYPE(ChebyshevPolynomial_t), INTENT(IN) :: poly
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(Matrix_ps), DIMENSION(:), ALLOCATABLE :: T_Powers
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Variables
    INTEGER :: degree
    INTEGER :: log2degree
    INTEGER :: counter

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Chebyshev Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Recursive")
       CALL WriteElement(key="Degree", int_value_in=degree)
       CALL PrintParameters(solver_parameters)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL CopyMatrix(InputMat,BalancedInput)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Construct The X Powers Array
    !! First, compute how many powers of two are necessary to compute.
    log2degree = 1
    DO WHILE (2**log2degree .LE. degree)
       log2degree = log2degree + 1
    END DO
    ALLOCATE(T_Powers(log2degree))

    !! Now compute those powers of two
    CALL CopyMatrix(Identity, T_Powers(1))
    IF (degree .EQ. 1) THEN
       CALL CopyMatrix(T_Powers(1), OutputMat)
    ELSE
       CALL CopyMatrix(BalancedInput,T_Powers(2))
       DO counter=3,log2degree
          CALL MatrixMultiply(T_Powers(counter-1), T_Powers(counter-1), &
               & T_Powers(counter), threshold_in=solver_parameters%threshold, &
               & alpha_in=REAL(2.0,NTREAL), memory_pool_in=pool)
          CALL IncrementMatrix(Identity, T_Powers(counter), &
               & alpha_in=REAL(-1.0,NTREAL))
       END DO

       !! Call Recursive
       CALL ComputeRecursive(T_Powers, poly, OutputMat, &
            &  pool, 1, solver_parameters)
    END IF
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DO counter=1,log2degree
       CALL DestructMatrix(T_Powers(counter))
    END DO
    DEALLOCATE(T_Powers)
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(BalancedInput)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE FactorizedChebyshevCompute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The workhorse routine for the factorized chebyshev computation function.
  !! @param[in] T_Powers the precomputed Chebyshev polynomials.
  !! @param[in] poly polynomial coefficients.
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[inout] pool the memory pool.
  !! @param[in] depth the depth of recursion.
  !! @param[in] solver_parameters parameters for the solver.
  RECURSIVE SUBROUTINE ComputeRecursive(T_Powers, poly, OutputMat, pool, &
       & depth, solver_parameters)
    !! Parameters
    TYPE(Matrix_ps), DIMENSION(:), INTENT(IN) :: T_Powers
    TYPE(ChebyshevPolynomial_t), INTENT(IN) :: poly
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    INTEGER, INTENT(in) :: depth
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !! Local Data
    INTEGER :: coefficient_midpoint
    INTEGER :: left_length, right_length
    INTEGER :: full_midpoint
    INTEGER :: counter
    TYPE(ChebyshevPolynomial_t) :: left_poly
    TYPE(ChebyshevPolynomial_t) :: right_poly
    TYPE(Matrix_ps) :: LeftMat
    TYPE(Matrix_ps) :: RightMat

    !! First Handle The Base Case
    IF (SIZE(poly%coefficients) .EQ. 2) THEN
       CALL CopyMatrix(T_Powers(1), OutputMat)
       CALL ScaleMatrix(OutputMat, poly%coefficients(1))
       CALL IncrementMatrix(T_Powers(2), OutputMat, &
            & alpha_in=poly%coefficients(2))
    ELSE
       !! Adjust the coefficients.
       coefficient_midpoint = SIZE(poly%coefficients)/2
       left_length = coefficient_midpoint
       right_length = SIZE(poly%coefficients) - coefficient_midpoint
       ALLOCATE(left_poly%coefficients(left_length))
       ALLOCATE(right_poly%coefficients(right_length))
       left_poly%coefficients(:) = poly%coefficients(:coefficient_midpoint)
       right_poly%coefficients(:) = poly%coefficients(coefficient_midpoint+1:)
       DO counter=2,SIZE(left_poly%coefficients)
          left_poly%coefficients(counter) = left_poly%coefficients(counter) - &
               & poly%coefficients(SIZE(poly%coefficients) - counter + 2)
       END DO

       !! Left recursion
       CALL ComputeRecursive(T_Powers, left_poly, LeftMat, pool, depth+1, &
            & solver_parameters)

       !! Right recursion
       full_midpoint = SIZE(T_Powers) - depth + 1
       CALL ComputeRecursive(T_Powers, right_poly, RightMat, pool, depth+1, &
            & solver_parameters)

       !! Sum Together
       CALL MatrixMultiply(T_Powers(full_midpoint), RightMat, &
            & OutputMat, threshold_in=solver_parameters%threshold, &
            & alpha_in=REAL(2.0,NTREAL), memory_pool_in=pool)

       CALL IncrementMatrix(LeftMat,OutputMat)
       CALL IncrementMatrix(T_Powers(full_midpoint), &
            & OutputMat, alpha_in=-1.0*right_poly%coefficients(1))

       !! Cleanup
       DEALLOCATE(left_poly%coefficients)
       DEALLOCATE(right_poly%coefficients)
       CALL DestructMatrix(LeftMat)
       CALL DestructMatrix(RightMat)
    END IF

  END SUBROUTINE ComputeRecursive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ChebyshevSolversModule
