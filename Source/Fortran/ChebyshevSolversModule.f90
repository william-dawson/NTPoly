!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing Matrix functions based on Chebyshev polynomials.
MODULE ChebyshevSolversModule
  USE DataTypesModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixModule
  USE FixedSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ProcessGridModule
  USE TimerModule
  USE mpi
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
    TYPE(ChebyshevPolynomial_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: degree

    ALLOCATE(this%coefficients(degree))
    this%coefficients = 0
  END SUBROUTINE ConstructChebyshevPolynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  !! @param[inout] this the polynomial to destruct.
  PURE SUBROUTINE DestructChebyshevPolynomial(this)
    !! Parameters
    TYPE(ChebyshevPolynomial_t), INTENT(inout) :: this
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
    TYPE(ChebyshevPolynomial_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: degree
    REAL(NTREAL), INTENT(in) :: coefficient

    this%coefficients(degree) = coefficient
  END SUBROUTINE SetChebyshevCoefficient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Chebyshev Polynomial of the matrix.
  !! This method uses the standard Chebyshev Polynomial expansion.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in] poly polynomial to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ChebyshevCompute(InputMat, OutputMat, poly, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: OutputMat
    TYPE(ChebyshevPolynomial_t), INTENT(in) :: poly
    TYPE(FixedSolverParameters_t), INTENT(in), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: BalancedInput
    TYPE(DistributedSparseMatrix_t) :: Tk
    TYPE(DistributedSparseMatrix_t) :: Tkminus1
    TYPE(DistributedSparseMatrix_t) :: Tkminus2
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Variables
    INTEGER :: degree
    INTEGER :: counter
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Chebyshev Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Standard")
       CALL WriteElement(key="Degree", int_value_in=degree)
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)
    CALL CopyDistributedSparseMatrix(InputMat,BalancedInput)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! First Term
    CALL CopyDistributedSparseMatrix(Identity,Tkminus2)
    IF (degree == 1) THEN
       CALL CopyDistributedSparseMatrix(Tkminus2,OutputMat)
       CALL ScaleDistributedSparseMatrix(OutputMat,poly%coefficients(1))
    ELSE
       CALL CopyDistributedSparseMatrix(BalancedInput,Tkminus1)
       CALL CopyDistributedSparseMatrix(Tkminus2,OutputMat)
       CALL ScaleDistributedSparseMatrix(OutputMat,poly%coefficients(1))
       CALL IncrementDistributedSparseMatrix(Tkminus1,OutputMat, &
            & alpha_in=poly%coefficients(2))
       IF (degree > 2) THEN
          CALL DistributedGemm(BalancedInput, Tkminus1, Tk, &
               & alpha_in=REAL(2.0,NTREAL), &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
          CALL IncrementDistributedSparseMatrix(Tkminus2,Tk,REAL(-1.0,NTREAL))
          CALL IncrementDistributedSparseMatrix(Tk, OutputMat, &
               & alpha_in=poly%coefficients(3))
          DO counter = 4, degree
             CALL CopyDistributedSparseMatrix(Tkminus1,Tkminus2)
             CALL CopyDistributedSparseMatrix(Tk,Tkminus1)
             CALL DistributedGemm(BalancedInput, Tkminus1, Tk, &
                  & alpha_in=REAL(2.0,NTREAL), &
                  & threshold_in=solver_parameters%threshold, &
                  & memory_pool_in=pool)
             CALL IncrementDistributedSparseMatrix(Tkminus2,Tk, &
                  & REAL(-1.0,NTREAL))
             CALL IncrementDistributedSparseMatrix(Tk, OutputMat, &
                  & alpha_in=poly%coefficients(counter))
          END DO
       END IF
    END IF
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(OutputMat,min_size,max_size)
       sparsity = GetSize(OutputMat) / &
            & (REAL(OutputMat%actual_matrix_dimension,KIND=NTREAL)**2)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
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
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(Tk)
    CALL DestructDistributedSparseMatrix(Tkminus1)
    CALL DestructDistributedSparseMatrix(Tkminus2)
    CALL DestructDistributedSparseMatrix(BalancedInput)
  END SUBROUTINE ChebyshevCompute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Chebyshev Polynomial of the matrix.
  !! This version first factors the Chebyshev Polynomial and computes the
  !! function using a divide and conquer algorithm. Based on a simplified
  !! version of the first method in \cite liang2003improved .
  !! @param[in]  InputMat the input matrix
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in]  poly polynomial to compute.
  !! @param[in]  solver_parameters_in parameters for the solver (optional).
  SUBROUTINE FactorizedChebyshevCompute(InputMat, OutputMat, poly, &
       & solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: InputMat
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: OutputMat
    TYPE(ChebyshevPolynomial_t), INTENT(IN) :: poly
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: BalancedInput
    TYPE(DistributedSparseMatrix_t), DIMENSION(:), ALLOCATABLE :: T_Powers
    TYPE(DistributedMatrixMemoryPool_t) :: pool
    !! Local Variables
    INTEGER :: degree
    INTEGER :: log2degree
    INTEGER :: counter
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF

    degree = SIZE(poly%coefficients)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Chebyshev Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Recursive")
       CALL WriteElement(key="Degree", int_value_in=degree)
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Initial values for matrices
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)
    CALL CopyDistributedSparseMatrix(InputMat,BalancedInput)

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
    CALL CopyDistributedSparseMatrix(Identity, T_Powers(1))
    IF (degree .EQ. 1) THEN
       CALL CopyDistributedSparseMatrix(T_Powers(1), OutputMat)
    ELSE
       CALL CopyDistributedSparseMatrix(BalancedInput,T_Powers(2))
       DO counter=3,log2degree
          CALL DistributedGemm(T_Powers(counter-1), T_Powers(counter-1), &
               & T_Powers(counter), threshold_in=solver_parameters%threshold, &
               & alpha_in=REAL(2.0,NTREAL), memory_pool_in=pool)
          CALL IncrementDistributedSparseMatrix(Identity, T_Powers(counter), &
               & alpha_in=REAL(-1.0,NTREAL))
       END DO

       !! Call Recursive
       CALL ComputeRecursive(T_Powers, poly, OutputMat, &
            &  pool, 1, solver_parameters)
    END IF
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(OutputMat,min_size,max_size)
       sparsity = GetSize(OutputMat) / &
            & (REAL(OutputMat%actual_matrix_dimension,KIND=NTREAL)**2)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
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
       CALL DestructDistributedSparseMatrix(T_Powers(counter))
    END DO
    DEALLOCATE(T_Powers)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(BalancedInput)
  END SUBROUTINE FactorizedChebyshevCompute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The workhorse routine for the factorized chebyshev computation function.
  !! @param[in] T_Powers the precomputed Chebyshev polynomials.
  !! @param[in] poly polynomial coefficients.
  !! @param[out] OutputMat = poly(InputMat)
  !! @param[in] depth_in the depth of recursion.
  !! @param[in] solver_parameters_in parameters for the solver.
  RECURSIVE SUBROUTINE ComputeRecursive(T_Powers, poly, OutputMat, pool, &
       & depth, solver_parameters)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), DIMENSION(:), INTENT(in) :: T_Powers
    TYPE(ChebyshevPolynomial_t), INTENT(in) :: poly
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: OutputMat
    INTEGER, INTENT(in) :: depth
    TYPE(FixedSolverParameters_t), INTENT(in) :: solver_parameters
    TYPE(DistributedMatrixMemoryPool_t), INTENT(inout) :: pool
    !! Local Data
    INTEGER :: coefficient_midpoint
    INTEGER :: left_length, right_length
    INTEGER :: full_midpoint
    INTEGER :: counter
    TYPE(ChebyshevPolynomial_t) :: left_poly
    TYPE(ChebyshevPolynomial_t) :: right_poly
    TYPE(DistributedSparseMatrix_t) :: LeftMat
    TYPE(DistributedSparseMatrix_t) :: RightMat

    !! First Handle The Base Case
    IF (SIZE(poly%coefficients) .EQ. 2) THEN
       CALL CopyDistributedSparseMatrix(T_Powers(1), OutputMat)
       CALL ScaleDistributedSparseMatrix(OutputMat, poly%coefficients(1))
       CALL IncrementDistributedSparseMatrix(T_Powers(2), OutputMat, &
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
       CALL DistributedGemm(T_Powers(full_midpoint), RightMat, &
            & OutputMat, threshold_in=solver_parameters%threshold, &
            & alpha_in=REAL(2.0,NTREAL), memory_pool_in=pool)

       CALL IncrementDistributedSparseMatrix(LeftMat,OutputMat)
       CALL IncrementDistributedSparseMatrix(T_Powers(full_midpoint), &
            & OutputMat, alpha_in=-1.0*right_poly%coefficients(1))

       !! Cleanup
       DEALLOCATE(left_poly%coefficients)
       DEALLOCATE(right_poly%coefficients)
       CALL DestructDistributedSparseMatrix(LeftMat)
       CALL DestructDistributedSparseMatrix(RightMat)
    END IF

  END SUBROUTINE ComputeRecursive
END MODULE ChebyshevSolversModule
