!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Matrix functions based on Chebyshev polynomials.
MODULE ChebyshevSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteHeader, EnterSubLog, ExitSubLog
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, FillMatrixIdentity, &
       & PrintMatrixInformation, ConstructEmptyMatrix, DestructMatrix, &
       & CopyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
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
  PUBLIC :: ConstructPolynomial
  PUBLIC :: DestructPolynomial
  PUBLIC :: SetCoefficient
  !! Solvers
  PUBLIC :: Compute
  PUBLIC :: FactorizedCompute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructPolynomial
     MODULE PROCEDURE ConstructPolynomial_cheby
  END INTERFACE ConstructPolynomial
  INTERFACE DestructPolynomial
     MODULE PROCEDURE DestructPolynomial_cheby
  END INTERFACE DestructPolynomial
  INTERFACE SetCoefficient
     MODULE PROCEDURE SetCoefficient_cheby
  END INTERFACE SetCoefficient
  INTERFACE Compute
     MODULE PROCEDURE Compute_cheby
  END INTERFACE Compute
  INTERFACE FactorizedCompute
     MODULE PROCEDURE FactorizedCompute_cheby
  END INTERFACE FactorizedCompute
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a Chebyshev polynomial object.
  PURE SUBROUTINE ConstructPolynomial_cheby(this, degree)
    !> The polynomial to construct.
    TYPE(ChebyshevPolynomial_t), INTENT(INOUT) :: this
    !> Degree of the polynomial.
    INTEGER, INTENT(IN) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0_NTREAL
  END SUBROUTINE ConstructPolynomial_cheby
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  PURE SUBROUTINE DestructPolynomial_cheby(this)
    !> The polynomial to destruct.
    TYPE(ChebyshevPolynomial_t), INTENT(INOUT) :: this
    IF (ALLOCATED(this%coefficients)) THEN
       DEALLOCATE(this%coefficients)
    END IF
  END SUBROUTINE DestructPolynomial_cheby
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set a coefficient of a Chebyshev polynomial.
  SUBROUTINE SetCoefficient_cheby(this, degree, coefficient)
    !> The polynomial to set.
    TYPE(ChebyshevPolynomial_t), INTENT(INOUT) :: this
    !> Degree for which to set the coefficient.
    INTEGER, INTENT(IN) :: degree
    !> Coefficient value.
    REAL(NTREAL), INTENT(IN) :: coefficient

    this%coefficients(degree) = coefficient
  END SUBROUTINE SetCoefficient_cheby
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Chebyshev Polynomial of the matrix.
  !> This method uses the standard Chebyshev Polynomial expansion.
  SUBROUTINE Compute_cheby(InputMat, OutputMat, poly, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> The Chebyshev polynomial to compute.
    TYPE(ChebyshevPolynomial_t), INTENT(IN) :: poly
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(Matrix_ps) :: Tk
    TYPE(Matrix_ps) :: Tkminus1
    TYPE(Matrix_ps) :: Tkminus2
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Variables
    INTEGER :: degree
    INTEGER :: II

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    degree = SIZE(poly%coefficients)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Chebyshev Solver")
       CALL EnterSubLog
       CALL WriteElement(key = "Method", VALUE = "Standard")
       CALL WriteElement(key = "Degree", VALUE = degree - 1)
       CALL PrintParameters(params)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL CopyMatrix(InputMat,BalancedInput)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in = pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & params%BalancePermutation, memorypool_in = pool)
    END IF

    !! First Term
    CALL CopyMatrix(Identity, Tkminus2)
    IF (degree == 1) THEN
       CALL CopyMatrix(Tkminus2, OutputMat)
       CALL ScaleMatrix(OutputMat, poly%coefficients(1))
    ELSE
       CALL CopyMatrix(BalancedInput, Tkminus1)
       CALL CopyMatrix(Tkminus2, OutputMat)
       CALL ScaleMatrix(OutputMat, poly%coefficients(1))
       CALL IncrementMatrix(Tkminus1, OutputMat, &
            & alpha_in = poly%coefficients(2))
       IF (degree .GT. 2) THEN
          CALL MatrixMultiply(BalancedInput, Tkminus1, Tk, &
               & alpha_in = 2.0_NTREAL, threshold_in = params%threshold, &
               & memory_pool_in = pool)
          CALL IncrementMatrix(Tkminus2, Tk, alpha_in = -1.0_NTREAL)
          CALL IncrementMatrix(Tk, OutputMat, &
               & alpha_in = poly%coefficients(3))
          DO II = 4, degree
             CALL CopyMatrix(Tkminus1, Tkminus2)
             CALL CopyMatrix(Tk, Tkminus1)
             CALL MatrixMultiply(BalancedInput, Tkminus1, Tk, &
                  & alpha_in = 2.0_NTREAL, threshold_in = params%threshold, &
                  & memory_pool_in = pool)
             CALL IncrementMatrix(Tkminus2, Tk, alpha_in = -1.0_NTREAL)
             CALL IncrementMatrix(Tk, OutputMat, &
                  & alpha_in = poly%coefficients(II))
          END DO
       END IF
    END IF
    IF (params%be_verbose) THEN
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & params%BalancePermutation, memorypool_in = pool)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(Tk)
    CALL DestructMatrix(Tkminus1)
    CALL DestructMatrix(Tkminus2)
    CALL DestructMatrix(BalancedInput)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE Compute_cheby
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Chebyshev Polynomial of the matrix.
  !> This version first factors the Chebyshev Polynomial and computes the
  !> function using a divide and conquer algorithm. Based on a simplified
  !> version of the first method in \cite liang2003improved .
  SUBROUTINE FactorizedCompute_cheby(InputMat, OutputMat, poly, &
       & solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> The Chebyshev polynomial to compute.
    TYPE(ChebyshevPolynomial_t), INTENT(IN) :: poly
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(Matrix_ps), DIMENSION(:), ALLOCATABLE :: T_Powers
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Variables
    INTEGER :: degree
    INTEGER :: log2degree
    INTEGER :: II

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    degree = SIZE(poly%coefficients)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Chebyshev Solver")
       CALL EnterSubLog
       CALL WriteElement(key = "Method", VALUE = "Recursive")
       CALL WriteElement(key = "Degree", VALUE = degree-1)
       CALL PrintParameters(params)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL CopyMatrix(InputMat, BalancedInput)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in = pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & params%BalancePermutation, memorypool_in = pool)
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
       CALL CopyMatrix(BalancedInput, T_Powers(2))
       DO II=3, log2degree
          CALL MatrixMultiply(T_Powers(II - 1), T_Powers(II - 1), &
               & T_Powers(II), threshold_in = params%threshold, &
               & alpha_in = 2.0_NTREAL, memory_pool_in = pool)
          CALL IncrementMatrix(Identity, T_Powers(II), alpha_in = -1.0_NTREAL)
       END DO
       !! Call Recursive
       CALL ComputeRecursive(T_Powers, poly, OutputMat, pool, 1, params)
    END IF
    IF (params%be_verbose) THEN
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & params%BalancePermutation, memorypool_in = pool)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DO II = 1, log2degree
       CALL DestructMatrix(T_Powers(II))
    END DO
    DEALLOCATE(T_Powers)
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(BalancedInput)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE FactorizedCompute_cheby
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The workhorse routine for the factorized chebyshev computation function.
  RECURSIVE SUBROUTINE ComputeRecursive(T_Powers, poly, OutputMat, pool, &
       & depth, params)
    !> The precomputed Chebyshev polynomials.
    TYPE(Matrix_ps), DIMENSION(:), INTENT(IN) :: T_Powers
    !> Polynomial coefficients.
    TYPE(ChebyshevPolynomial_t), INTENT(IN) :: poly
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> The depth of recursion.
    INTEGER, INTENT(in) :: depth
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !> The memory pool.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !! Local Data
    INTEGER :: coefficient_midpoint
    INTEGER :: left_length, right_length
    INTEGER :: full_midpoint
    TYPE(ChebyshevPolynomial_t) :: left_poly
    TYPE(ChebyshevPolynomial_t) :: right_poly
    TYPE(Matrix_ps) :: LeftMat
    TYPE(Matrix_ps) :: RightMat
    INTEGER :: II

    !! First Handle The Base Case
    IF (SIZE(poly%coefficients) .EQ. 1) THEN
       CALL CopyMatrix(T_Powers(1), OutputMat)
       CALL ScaleMatrix(OutputMat, poly%coefficients(1))
    ELSE IF (SIZE(poly%coefficients) .EQ. 2) THEN
       CALL CopyMatrix(T_Powers(1), OutputMat)
       CALL ScaleMatrix(OutputMat, poly%coefficients(1))
       CALL IncrementMatrix(T_Powers(2), OutputMat, &
            & alpha_in=poly%coefficients(2))
    ELSE
       !! Adjust the coefficients.
       coefficient_midpoint = SIZE(poly%coefficients) / 2
       left_length = coefficient_midpoint
       right_length = SIZE(poly%coefficients) - coefficient_midpoint
       ALLOCATE(left_poly%coefficients(left_length))
       ALLOCATE(right_poly%coefficients(right_length))
       left_poly%coefficients(:) = poly%coefficients(:coefficient_midpoint)
       right_poly%coefficients(:) = poly%coefficients(coefficient_midpoint+1:)
       DO II = 2, SIZE(left_poly%coefficients)
          left_poly%coefficients(II) = left_poly%coefficients(II) - &
               & poly%coefficients(SIZE(poly%coefficients) - II + 2)
       END DO

       !! Left recursion
       CALL ComputeRecursive(T_Powers, left_poly, LeftMat, pool, depth + 1, &
            & params)

       !! Right recursion
       full_midpoint = SIZE(T_Powers) - depth + 1
       CALL ComputeRecursive(T_Powers, right_poly, RightMat, pool, depth + 1, &
            & params)

       !! Sum Together
       CALL MatrixMultiply(T_Powers(full_midpoint), RightMat, &
            & OutputMat, threshold_in = params%threshold, &
            & alpha_in = 2.0_NTREAL, memory_pool_in = pool)

       CALL IncrementMatrix(LeftMat, OutputMat)
       CALL IncrementMatrix(T_Powers(full_midpoint), &
            & OutputMat, alpha_in = -1.0*right_poly%coefficients(1))

       !! Cleanup
       DEALLOCATE(left_poly%coefficients)
       DEALLOCATE(right_poly%coefficients)
       CALL DestructMatrix(LeftMat)
       CALL DestructMatrix(RightMat)
    END IF

  END SUBROUTINE ComputeRecursive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ChebyshevSolversModule
