!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing matrix functions based on Hermite polynomials.
!! The Physicist variety.
MODULE HermiteSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, WriteHeader
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix, MatrixMultiply, &
       & ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, FillMatrixIdentity, PrintMatrixInformation
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype that represents a Hermite polynomial.
  TYPE, PUBLIC :: HermitePolynomial_t
     !> Coefficients of the polynomial.
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: coefficients
  END TYPE HermitePolynomial_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Polynomial type
  PUBLIC :: ConstructPolynomial
  PUBLIC :: DestructPolynomial
  PUBLIC :: SetCoefficient
  !! Solvers
  PUBLIC :: Compute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructPolynomial
     MODULE PROCEDURE ConstructPolynomial_horner
  END INTERFACE ConstructPolynomial
  INTERFACE DestructPolynomial
     MODULE PROCEDURE DestructPolynomial_horner
  END INTERFACE DestructPolynomial
  INTERFACE SetCoefficient
     MODULE PROCEDURE SetCoefficient_horner
  END INTERFACE SetCoefficient
  INTERFACE Compute
     MODULE PROCEDURE Compute_horner
  END INTERFACE Compute
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a Hermite polynomial object.
  PURE SUBROUTINE ConstructPolynomial_horner(this, degree)
    !> The polynomial to construct.
    TYPE(HermitePolynomial_t), INTENT(inout) :: this
    !> The degree of the polynomial.
    INTEGER, INTENT(in) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0.0_NTREAL
  END SUBROUTINE ConstructPolynomial_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Hermite polynomial object.
  PURE SUBROUTINE DestructPolynomial_horner(this)
    !> The polynomial to destruct.
    TYPE(HermitePolynomial_t), INTENT(inout) :: this

    IF (ALLOCATED(this%coefficients)) THEN
       DEALLOCATE(this%coefficients)
    END IF
  END SUBROUTINE DestructPolynomial_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set a coefficient of a Hermite polynomial.
  SUBROUTINE SetCoefficient_horner(this, degree, coefficient)
    !> The polynomial to set.
    TYPE(HermitePolynomial_t), INTENT(inout) :: this
    !> The degree for which to set the coefficient.
    INTEGER, INTENT(in) :: degree
    !> Coefficient value to set.
    REAL(NTREAL), INTENT(in) :: coefficient

    this%coefficients(degree) = coefficient
  END SUBROUTINE SetCoefficient_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Hermite Polynomial of the matrix.
  !> This method uses the standard Hermite Polynomial expansion.
  SUBROUTINE Compute_horner(InputMat, OutputMat, poly, solver_parameters_in)
    !> The input matrix.
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    !> Polynomial to compute.
    TYPE(HermitePolynomial_t), INTENT(in) :: poly
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(Matrix_ps) :: Hk
    TYPE(Matrix_ps) :: Hkminus1
    TYPE(Matrix_ps) :: Hkplus1
    TYPE(Matrix_ps) :: Hkprime
    TYPE(MatrixMemoryPool_p) :: pool
    !! Local Variables
    INTEGER :: degree
    INTEGER :: II

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Hermite Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Standard")
       CALL WriteElement(key="Degree", VALUE=degree-1)
       CALL PrintParameters(params)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL CopyMatrix(InputMat, BalancedInput)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Recursive expansion
    CALL CopyMatrix(Identity, Hkminus1)
    CALL CopyMatrix(Hkminus1, OutputMat)
    CALL ScaleMatrix(OutputMat, poly%coefficients(1))
    IF (degree .GT. 1) THEN
       CALL CopyMatrix(BalancedInput, Hk)
       CALL ScaleMatrix(Hk, REAL(2.0,KIND=NTREAL))
       CALL IncrementMatrix(Hk, OutputMat, &
            & alpha_in=poly%coefficients(2))
       IF (degree .GT. 2) THEN
          CALL CopyMatrix(Hkminus1, Hkprime)
          CALL ScaleMatrix(Hkprime, REAL(2.0,NTREAL))
          DO II = 3, degree
             CALL MatrixMultiply(BalancedInput, Hk, Hkplus1, &
                  & alpha_in=REAL(2.0,NTREAL), &
                  & threshold_in=params%threshold, &
                  & memory_pool_in=pool)
             CALL IncrementMatrix(Hkprime, Hkplus1, &
                  & alpha_in=REAL(-1.0,NTREAL))
             CALL CopyMatrix(Hk, Hkprime)
             CALL ScaleMatrix(Hkprime, &
                  & REAL(2*(II-1),KIND=NTREAL))
             CALL CopyMatrix(Hk, Hkminus1)
             CALL CopyMatrix(Hkplus1, Hk)
             CALL IncrementMatrix(Hk, OutputMat, &
                  & alpha_in=poly%coefficients(II))
          END DO
       END IF
    END IF
    IF (params%be_verbose) THEN
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(Hk)
    CALL DestructMatrix(Hkminus1)
    CALL DestructMatrix(Hkplus1)
    CALL DestructMatrix(Hkprime)
    CALL DestructMatrix(BalancedInput)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE Compute_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE HermiteSolversModule
