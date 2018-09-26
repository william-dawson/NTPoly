!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing General Matrix Polynomials.
MODULE PolynomialSolversModule
  USE DataTypesModule
  USE LoadBalancerModule
  USE LoggingModule
  USE PMatrixMemoryPoolModule
  USE PSMatrixAlgebraModule
  USE PSMatrixModule
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE TimerModule
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
  END INTERFACE
  INTERFACE DestructPolynomial
     MODULE PROCEDURE DestructPolynomial_stand
  END INTERFACE
  INTERFACE SetCoefficient
     MODULE PROCEDURE SetCoefficient_stand
  END INTERFACE
  INTERFACE Compute
     MODULE PROCEDURE Compute_stand
  END INTERFACE
  INTERFACE FactorizedCompute
     MODULE PROCEDURE FactorizedCompute_stand
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a polynomial.
  PURE SUBROUTINE ConstructPolynomial_stand(this, degree)
    !> The polynomial to construct.
    TYPE(Polynomial_t), INTENT(INOUT) :: this
    !> The degree of the polynomial.
    INTEGER, INTENT(IN) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0
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
  !> Compute A Matrix Polynomial Using Horner's Method.
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
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(Matrix_ps) :: Temporary
    INTEGER :: degree
    INTEGER :: counter
    TYPE(MatrixMemoryPool_p) :: pool

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Polynomial Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Horner")
       CALL PrintParameters(solver_parameters)
       CALL WriteElement(key="Degree", int_value_in=degree-1)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL ConstructEmptyMatrix(Temporary, InputMat)
    CALL CopyMatrix(InputMat,BalancedInput)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF
    CALL CopyMatrix(Identity,OutputMat)

    IF (SIZE(poly%coefficients) .EQ. 1) THEN
       CALL ScaleMatrix(OutputMat, poly%coefficients(degree))
    ELSE
       CALL ScaleMatrix(OutputMat,poly%coefficients(degree-1))
       CALL IncrementMatrix(BalancedInput,OutputMat, &
            & poly%coefficients(degree))
       DO counter = degree-2,1,-1
          CALL MatrixMultiply(BalancedInput,OutputMat,Temporary, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
          CALL CopyMatrix(Temporary,OutputMat)
          CALL IncrementMatrix(Identity, &
               & OutputMat, alpha_in=poly%coefficients(counter))
       END DO
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
    CALL DestructMatrix(Temporary)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE Compute_stand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using Paterson and Stockmeyer's method.
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
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps), DIMENSION(:), ALLOCATABLE :: x_powers
    TYPE(Matrix_ps) :: Bk
    TYPE(Matrix_ps) :: Xs
    TYPE(Matrix_ps) :: Temp
    INTEGER :: degree
    INTEGER :: m_value, s_value, r_value
    INTEGER :: k_value
    INTEGER :: counter
    INTEGER :: c_index
    TYPE(MatrixMemoryPool_p) :: pool

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    !! Parameters for splitting up polynomial.
    degree = SIZE(poly%coefficients)
    m_value = degree-1
    s_value = INT(SQRT(REAL(m_value)))
    r_value = m_value/s_value

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Polynomial Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Paterson Stockmeyer")
       CALL WriteCitation("paterson1973number")
       CALL PrintParameters(solver_parameters)
       CALL WriteElement(key="Degree", int_value_in=degree-1)
    END IF

    ALLOCATE(x_powers(s_value+1))

    !! Initial values for matrices
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)

    !! Create the X Powers
    CALL ConstructEmptyMatrix(x_powers(1), InputMat)
    CALL FillMatrixIdentity(x_powers(1))
    DO counter=1,s_value+1-1
       CALL MatrixMultiply(InputMat,x_powers(counter-1+1),x_powers(counter+1),&
            & memory_pool_in=pool)
    END DO
    CALL CopyMatrix(x_powers(s_value+1),Xs)

    !! S_k = bmX
    CALL CopyMatrix(Identity,Bk)
    CALL ScaleMatrix(Bk, poly%coefficients(s_value*r_value+1))
    DO counter=1,m_value-s_value*r_value+1-1
       c_index = s_value*r_value + counter
       CALL IncrementMatrix(x_powers(counter+1),Bk, &
            & alpha_in=poly%coefficients(c_index+1))
    END DO
    CALL MatrixMultiply(Bk,Xs,OutputMat, memory_pool_in=pool)

    !! S_k += bmx + bm-1I
    k_value = r_value - 1
    CALL CopyMatrix(Identity,Bk)
    CALL ScaleMatrix(Bk,poly%coefficients(s_value*k_value+1))
    DO counter=1,s_value-1+1-1
       c_index = s_value*k_value + counter
       CALL IncrementMatrix(x_powers(counter+1),Bk, &
            & alpha_in=poly%coefficients(c_index+1))
    END DO
    CALL IncrementMatrix(Bk,OutputMat)

    !! Loop over the rest.
    DO k_value=r_value-2,-1+1,-1
       CALL CopyMatrix(Identity,Bk)
       CALL ScaleMatrix(Bk, &
            & poly%coefficients(s_value*k_value+1))
       DO counter=1,s_value-1+1-1
          c_index = s_value*k_value + counter
          CALL IncrementMatrix(x_powers(counter+1),Bk, &
               & alpha_in=poly%coefficients(c_index+1))
       END DO
       CALL MatrixMultiply(Xs,OutputMat,Temp)
       CALL CopyMatrix(Temp,OutputMat)
       CALL IncrementMatrix(Bk,OutputMat)
    END DO

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DO counter=1,s_value+1
       CALL DestructMatrix(x_powers(counter))
    END DO
    DEALLOCATE(x_powers)
    CALL DestructMatrix(Bk)
    CALL DestructMatrix(Xs)
    CALL DestructMatrix(Temp)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE FactorizedCompute_stand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PolynomialSolversModule
