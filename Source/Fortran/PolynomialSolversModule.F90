!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing General Matrix Polynomials.
MODULE PolynomialSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix, MatrixMultiply, &
       & ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, DestructMatrix, FillMatrixIdentity, &
       & ConstructEmptyMatrix, CopyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype that represents a polynomial.
  TYPE, PUBLIC :: Polynomial_t
     !> Coefficients of the polynomial.
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: coefficients
  END TYPE Polynomial_t
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
     MODULE PROCEDURE ConstructPolynomial_stand
  END INTERFACE ConstructPolynomial
  INTERFACE DestructPolynomial
     MODULE PROCEDURE DestructPolynomial_stand
  END INTERFACE DestructPolynomial
  INTERFACE SetCoefficient
     MODULE PROCEDURE SetCoefficient_stand
  END INTERFACE SetCoefficient
  INTERFACE Compute
     MODULE PROCEDURE Compute_stand
  END INTERFACE Compute
  INTERFACE FactorizedCompute
     MODULE PROCEDURE FactorizedCompute_stand
  END INTERFACE FactorizedCompute
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a polynomial.
  PURE SUBROUTINE ConstructPolynomial_stand(this, degree)
    !> The polynomial to construct.
    TYPE(Polynomial_t), INTENT(INOUT) :: this
    !> The degree of the polynomial.
    INTEGER, INTENT(IN) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0_NTREAL
  END SUBROUTINE ConstructPolynomial_stand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  PURE SUBROUTINE DestructPolynomial_stand(this)
    !> The polynomial to destruct.
    TYPE(Polynomial_t), INTENT(INOUT) :: this
    IF (ALLOCATED(this%coefficients)) THEN
       DEALLOCATE(this%coefficients)
    END IF
  END SUBROUTINE DestructPolynomial_stand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set coefficient of a polynomial.
  SUBROUTINE SetCoefficient_stand(this, degree, coefficient)
    !> The polynomial to set.
    TYPE(Polynomial_t), INTENT(INOUT) :: this
    !> Degree for which to set the coefficient.
    INTEGER, INTENT(IN) :: degree
    !> Coefficient value.
    REAL(NTREAL), INTENT(IN) :: coefficient

    this%coefficients(degree) = coefficient
  END SUBROUTINE SetCoefficient_stand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using the method of Horner.
  SUBROUTINE Compute_stand(InputMat, OutputMat, poly, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Polynomial to compute.
    TYPE(Polynomial_t), INTENT(IN) :: poly
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(Matrix_ps) :: Temporary
    INTEGER :: degree
    INTEGER :: II
    TYPE(MatrixMemoryPool_p) :: pool

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Polynomial Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Horner")
       CALL PrintParameters(params)
       CALL WriteElement(key="Degree", VALUE=degree-1)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL ConstructEmptyMatrix(Temporary, InputMat)
    CALL CopyMatrix(InputMat, BalancedInput)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF
    CALL CopyMatrix(Identity, OutputMat)

    IF (SIZE(poly%coefficients) .EQ. 1) THEN
       CALL ScaleMatrix(OutputMat, poly%coefficients(degree))
    ELSE
       CALL ScaleMatrix(OutputMat,poly%coefficients(degree-1))
       CALL IncrementMatrix(BalancedInput,OutputMat, &
            & poly%coefficients(degree))
       DO II = degree-2, 1, -1
          CALL MatrixMultiply(BalancedInput,OutputMat,Temporary, &
               & threshold_in=params%threshold, memory_pool_in=pool)
          CALL CopyMatrix(Temporary,OutputMat)
          CALL IncrementMatrix(Identity, &
               & OutputMat, alpha_in=poly%coefficients(II))
       END DO
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
    CALL DestructMatrix(Temporary)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE Compute_stand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using The Paterson and Stockmeyer method.
  !> This method first factors the polynomial to reduce the number of
  !> matrix multiplies required.
  SUBROUTINE FactorizedCompute_stand(InputMat, OutputMat, poly, &
       & solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> The polynomial to compute.
    TYPE(Polynomial_t), INTENT(IN) :: poly
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps), DIMENSION(:), ALLOCATABLE :: x_powers
    TYPE(Matrix_ps) :: Bk
    TYPE(Matrix_ps) :: Xs
    TYPE(Matrix_ps) :: Temp
    INTEGER :: degree
    INTEGER :: m_value, s_value, r_value
    INTEGER :: k_value
    INTEGER :: II
    INTEGER :: c_index
    TYPE(MatrixMemoryPool_p) :: pool

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    !! Parameters for splitting up polynomial.
    degree = SIZE(poly%coefficients)
    m_value = degree-1
    s_value = INT(SQRT(REAL(m_value)))
    r_value = m_value/s_value

    IF (params%be_verbose) THEN
       CALL WriteHeader("Polynomial Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Paterson Stockmeyer")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("paterson1973number")
       CALL ExitSubLog
       CALL PrintParameters(params)
       CALL WriteElement(key="Degree", VALUE=degree-1)
    END IF

    ALLOCATE(x_powers(s_value+1))

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)

    !! Create the X Powers
    CALL ConstructEmptyMatrix(x_powers(1), InputMat)
    CALL FillMatrixIdentity(x_powers(1))
    DO II = 1, s_value+1-1
       CALL MatrixMultiply(InputMat,x_powers(II-1+1),x_powers(II+1),&
            & memory_pool_in=pool)
    END DO
    CALL CopyMatrix(x_powers(s_value+1),Xs)

    !! S_k = bmX
    CALL CopyMatrix(Identity,Bk)
    CALL ScaleMatrix(Bk, poly%coefficients(s_value*r_value+1))
    DO II = 1, m_value-s_value*r_value+1-1
       c_index = s_value*r_value + II
       CALL IncrementMatrix(x_powers(II+1),Bk, &
            & alpha_in=poly%coefficients(c_index+1))
    END DO
    CALL MatrixMultiply(Bk, Xs, OutputMat, memory_pool_in=pool)

    !! S_k += bmx + bm-1I
    k_value = r_value - 1
    CALL CopyMatrix(Identity, Bk)
    CALL ScaleMatrix(Bk, poly%coefficients(s_value*k_value+1))
    DO II = 1, s_value-1+1-1
       c_index = s_value*k_value + II
       CALL IncrementMatrix(x_powers(II+1),Bk, &
            & alpha_in=poly%coefficients(c_index+1))
    END DO
    CALL IncrementMatrix(Bk,OutputMat)

    !! Loop over the rest.
    DO k_value=r_value-2,-1+1,-1
       CALL CopyMatrix(Identity,Bk)
       CALL ScaleMatrix(Bk, poly%coefficients(s_value*k_value+1))
       DO II=1,s_value-1+1-1
          c_index = s_value*k_value + II
          CALL IncrementMatrix(x_powers(II+1),Bk, &
               & alpha_in=poly%coefficients(c_index+1))
       END DO
       CALL MatrixMultiply(Xs,OutputMat,Temp)
       CALL CopyMatrix(Temp,OutputMat)
       CALL IncrementMatrix(Bk,OutputMat)
    END DO

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DO II = 1, s_value+1
       CALL DestructMatrix(x_powers(II))
    END DO
    DEALLOCATE(x_powers)
    CALL DestructMatrix(Bk)
    CALL DestructMatrix(Xs)
    CALL DestructMatrix(Temp)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE FactorizedCompute_stand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PolynomialSolversModule
