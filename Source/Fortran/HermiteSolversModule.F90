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
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE MPI
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
  END INTERFACE
  INTERFACE DestructPolynomial
     MODULE PROCEDURE DestructPolynomial_horner
  END INTERFACE
  INTERFACE SetCoefficient
     MODULE PROCEDURE SetCoefficient_horner
  END INTERFACE
  INTERFACE Compute
     MODULE PROCEDURE Compute_horner
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a Hermite polynomial object.
  !! @param[inout] this the polynomial to construct.
  !! @param[in] degree of the polynomial.
  PURE SUBROUTINE ConstructPolynomial_horner(this, degree)
    !! Parameters
    TYPE(HermitePolynomial_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0
  END SUBROUTINE ConstructPolynomial_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Hermite polynomial object.
  !! @param[inout] this the polynomial to destruct.
  PURE SUBROUTINE DestructPolynomial_horner(this)
    !! Parameters
    TYPE(HermitePolynomial_t), INTENT(inout) :: this
    IF (ALLOCATED(this%coefficients)) THEN
       DEALLOCATE(this%coefficients)
    END IF
  END SUBROUTINE DestructPolynomial_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set a coefficient of a Hermite polynomial.
  !! @param[inout] this the polynomial to set.
  !! @param[in] degree for which to set the coefficient.
  !! @param[in] coefficient value.
  SUBROUTINE SetCoefficient_horner(this, degree, coefficient)
    !! Parameters
    TYPE(HermitePolynomial_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: degree
    REAL(NTREAL), INTENT(in) :: coefficient

    this%coefficients(degree) = coefficient
  END SUBROUTINE SetCoefficient_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Hermite Polynomial of the matrix.
  !! This method uses the standard Hermite Polynomial expansion.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in] poly polynomial to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE Compute_horner(InputMat, OutputMat, poly, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(HermitePolynomial_t), INTENT(in) :: poly
    TYPE(SolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
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
    INTEGER :: counter

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Hermite Solver")
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
          DO counter = 3, degree
             CALL MatrixMultiply(BalancedInput, Hk, Hkplus1, &
                  & alpha_in=REAL(2.0,NTREAL), &
                  & threshold_in=solver_parameters%threshold, &
                  & memory_pool_in=pool)
             CALL IncrementMatrix(Hkprime, Hkplus1, &
                  & alpha_in=REAL(-1.0,NTREAL))
             CALL CopyMatrix(Hk, Hkprime)
             CALL ScaleMatrix(Hkprime, &
                  & REAL(2*(counter-1),KIND=NTREAL))
             CALL CopyMatrix(Hk, Hkminus1)
             CALL CopyMatrix(Hkplus1, Hk)
             CALL IncrementMatrix(Hk, OutputMat, &
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
    CALL DestructMatrix(Hk)
    CALL DestructMatrix(Hkminus1)
    CALL DestructMatrix(Hkplus1)
    CALL DestructMatrix(Hkprime)
    CALL DestructMatrix(BalancedInput)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE Compute_horner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE HermiteSolversModule
