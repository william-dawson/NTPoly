!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing General Matrix Polynomials.
MODULE PolynomialSolversModule
  USE DataTypesModule
  USE MatrixMemoryPoolDModule
  USE MatrixDSAlgebraModule
  USE MatrixDSModule
  USE FixedSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE TimerModule
  USE MPI
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
  PUBLIC :: HornerCompute
  PUBLIC :: PatersonStockmeyerCompute
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a polynomial.
  !! @param[inout] this the polynomial to construct.
  !! @param[in] degree of the polynomial.
  PURE SUBROUTINE ConstructPolynomial(this, degree)
    !! Parameters
    TYPE(Polynomial_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0
  END SUBROUTINE ConstructPolynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  !! @param[inout] this the polynomial to destruct.
  PURE SUBROUTINE DestructPolynomial(this)
    !! Parameters
    TYPE(Polynomial_t), INTENT(inout) :: this
    IF (ALLOCATED(this%coefficients)) THEN
       DEALLOCATE(this%coefficients)
    END IF
  END SUBROUTINE DestructPolynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set coefficient of a polynomial.
  !! @param[inout] this the polynomial to set.
  !! @param[in] degree for which to set the coefficient.
  !! @param[in] coefficient value.
  SUBROUTINE SetCoefficient(this, degree, coefficient)
    !! Parameters
    TYPE(Polynomial_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: degree
    REAL(NTREAL), INTENT(in) :: coefficient

    this%coefficients(degree) = coefficient
  END SUBROUTINE SetCoefficient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using Horner's Method.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in] poly polynomial to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE HornerCompute(InputMat, OutputMat, poly, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ds), INTENT(in)  :: InputMat
    TYPE(Matrix_ds), INTENT(inout) :: OutputMat
    TYPE(Polynomial_t), INTENT(in) :: poly
    TYPE(FixedSolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ds) :: Identity
    TYPE(Matrix_ds) :: BalancedInput
    TYPE(Matrix_ds) :: Temporary
    INTEGER :: degree
    INTEGER :: counter
    TYPE(MatrixMemoryPool_d) :: pool

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Polynomial Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Horner")
       CALL PrintFixedSolverParameters(solver_parameters)
       CALL WriteElement(key="Degree", int_value_in=degree)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyMatrixDS(Identity, &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL FillDistributedIdentity(Identity)
    CALL ConstructEmptyMatrixDS(Temporary, &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL CopyMatrixDS(InputMat,BalancedInput)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF
    CALL CopyMatrixDS(Identity,OutputMat)

    IF (SIZE(poly%coefficients) .EQ. 1) THEN
       CALL ScaleDistributedSparseMatrix(OutputMat, poly%coefficients(degree))
    ELSE
       CALL ScaleDistributedSparseMatrix(OutputMat,poly%coefficients(degree-1))
       CALL IncrementDistributedSparseMatrix(BalancedInput,OutputMat, &
            & poly%coefficients(degree))
       DO counter = degree-2,1,-1
          CALL DistributedGemm(BalancedInput,OutputMat,Temporary, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
          CALL CopyMatrixDS(Temporary,OutputMat)
          CALL IncrementDistributedSparseMatrix(Identity, &
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
    CALL DestructDistributedSparseMatrix(Temporary)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructMatrixMemoryPoolD(pool)
  END SUBROUTINE HornerCompute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using Paterson and Stockmeyer's method.
  !! This method first factors the polynomial to reduce the number of
  !! matrix multiplies required.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in] poly polynomial to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE PatersonStockmeyerCompute(InputMat, OutputMat, poly, &
       & solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ds), INTENT(in)  :: InputMat
    TYPE(Matrix_ds), INTENT(inout) :: OutputMat
    TYPE(Polynomial_t), INTENT(in) :: poly
    TYPE(FixedSolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ds) :: Identity
    TYPE(Matrix_ds), DIMENSION(:), ALLOCATABLE :: x_powers
    TYPE(Matrix_ds) :: Bk
    TYPE(Matrix_ds) :: Xs
    TYPE(Matrix_ds) :: Temp
    INTEGER :: degree
    INTEGER :: m_value, s_value, r_value
    INTEGER :: k_value
    INTEGER :: counter
    INTEGER :: c_index
    TYPE(MatrixMemoryPool_d) :: pool

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
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
       CALL PrintFixedSolverParameters(solver_parameters)
       CALL WriteElement(key="Degree", int_value_in=degree)
    END IF

    ALLOCATE(x_powers(s_value+1))

    !! Initial values for matrices
    CALL ConstructEmptyMatrixDS(Identity, &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL FillDistributedIdentity(Identity)

    !! Create the X Powers
    CALL ConstructEmptyMatrixDS(x_powers(1), &
         & InputMat%actual_matrix_dimension, InputMat%process_grid)
    CALL FillDistributedIdentity(x_powers(1))
    DO counter=1,s_value+1-1
       CALL DistributedGemm(InputMat,x_powers(counter-1+1),x_powers(counter+1),&
            & memory_pool_in=pool)
    END DO
    CALL CopyMatrixDS(x_powers(s_value+1),Xs)

    !! S_k = bmX
    CALL CopyMatrixDS(Identity,Bk)
    CALL ScaleDistributedSparseMatrix(Bk, poly%coefficients(s_value*r_value+1))
    DO counter=1,m_value-s_value*r_value+1-1
       c_index = s_value*r_value + counter
       CALL IncrementDistributedSparseMatrix(x_powers(counter+1),Bk, &
            & alpha_in=poly%coefficients(c_index+1))
    END DO
    CALL DistributedGemm(Bk,Xs,OutputMat, memory_pool_in=pool)

    !! S_k += bmx + bm-1I
    k_value = r_value - 1
    CALL CopyMatrixDS(Identity,Bk)
    CALL ScaleDistributedSparseMatrix(Bk,poly%coefficients(s_value*k_value+1))
    DO counter=1,s_value-1+1-1
       c_index = s_value*k_value + counter
       CALL IncrementDistributedSparseMatrix(x_powers(counter+1),Bk, &
            & alpha_in=poly%coefficients(c_index+1))
    END DO
    CALL IncrementDistributedSparseMatrix(Bk,OutputMat)

    !! Loop over the rest.
    DO k_value=r_value-2,-1+1,-1
       CALL CopyMatrixDS(Identity,Bk)
       CALL ScaleDistributedSparseMatrix(Bk, &
            & poly%coefficients(s_value*k_value+1))
       DO counter=1,s_value-1+1-1
          c_index = s_value*k_value + counter
          CALL IncrementDistributedSparseMatrix(x_powers(counter+1),Bk, &
               & alpha_in=poly%coefficients(c_index+1))
       END DO
       CALL DistributedGemm(Xs,OutputMat,Temp)
       CALL CopyMatrixDS(Temp,OutputMat)
       CALL IncrementDistributedSparseMatrix(Bk,OutputMat)
    END DO

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DO counter=1,s_value+1
       CALL DestructDistributedSparseMatrix(x_powers(counter))
    END DO
    DEALLOCATE(x_powers)
    CALL DestructDistributedSparseMatrix(Bk)
    CALL DestructDistributedSparseMatrix(Xs)
    CALL DestructDistributedSparseMatrix(Temp)
    CALL DestructMatrixMemoryPoolD(pool)
  END SUBROUTINE PatersonStockmeyerCompute
END MODULE PolynomialSolversModule
