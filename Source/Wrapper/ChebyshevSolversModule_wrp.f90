!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the Chebyshev Solvers Module
MODULE ChebyshevSolversModule_wrp
  USE ChebyshevSolversModule, ONLY : ChebyshevPolynomial_t, &
       & ConstructChebyshevPolynomial, DestructChebyshevPolynomial, &
       & SetChebyshevCoefficient, ChebyshevCompute, FactorizedChebyshevCompute
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE FixedSolversModule_wrp, ONLY : FixedSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the polynomial data type.
  TYPE, PUBLIC :: ChebyshevPolynomial_wrp
     !> Actual data.
     TYPE(ChebyshevPolynomial_t), POINTER :: DATA
  END TYPE ChebyshevPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Chebyshev Polynomial type.
  PUBLIC :: ConstructChebyshevPolynomial_wrp
  PUBLIC :: DestructChebyshevPolynomial_wrp
  PUBLIC :: SetChebyshevCoefficient_wrp
  !! Solvers.
  PUBLIC :: ChebyshevCompute_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the empty polynomial constructor.
  !! @param[out] ih_this handle to the polynomial being created.
  !! @param[in] degree the degree of the polynomial.
  PURE SUBROUTINE ConstructChebyshevPolynomial_wrp(ih_this, degree) &
       & bind(c,name="ConstructChebyshevPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: degree
    TYPE(ChebyshevPolynomial_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructChebyshevPolynomial(h_this%data,degree)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructChebyshevPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  !! @param[inout] ih_this handle to the polynomial to free up.
  PURE SUBROUTINE DestructChebyshevPolynomial_wrp(ih_this) &
       & bind(c,name="DestructChebyshevPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(ChebyshevPolynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructChebyshevPolynomial(h_this%data)
    DEALLOCATE(h_this%data)
  END SUBROUTINE DestructChebyshevPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set coefficient of a polynomial.
  !! @param[inout] ih_this handle to the polynomial to set.
  !! @param[in] degree for which to set the coefficient.
  !! @param[in] coefficient value.
  SUBROUTINE SetChebyshevCoefficient_wrp(ih_this, degree, coefficient) &
       & bind(c,name="SetChebyshevCoefficient_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: degree
    REAL(NTREAL), INTENT(in) :: coefficient
    TYPE(ChebyshevPolynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetChebyshevCoefficient(h_this%data, degree, coefficient)
  END SUBROUTINE SetChebyshevCoefficient_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Chebyshev Polynomial.
  !! @param[in] ih_InputMat the input matrix
  !! @param[out] ih_OutputMat = poly(InputMat)
  !! @param[in] ih_polynomial polynomial to compute.
  !! @param[in] ih_solver_parameters parameters for the solver.
  SUBROUTINE ChebyshevCompute_wrp(ih_InputMat, ih_OutputMat, ih_polynomial, &
       & ih_solver_parameters) bind(c,name="ChebyshevCompute_wrp")
    INTEGER(kind=c_int), INTENT(in)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(ChebyshevPolynomial_wrp)     :: h_polynomial
    TYPE(FixedSolverParameters_wrp)   :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_polynomial = TRANSFER(ih_polynomial, h_polynomial)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ChebyshevCompute(h_InputMat%data, h_OutputMat%data, &
         & h_polynomial%data, h_solver_parameters%data)
  END SUBROUTINE ChebyshevCompute_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Chebyshev Polynomial By Factorization.
  !! @param[in] ih_InputMat the input matrix
  !! @param[out] ih_OutputMat = poly(InputMat)
  !! @param[in] ih_polynomial polynomial to compute.
  !! @param[in] ih_solver_parameters parameters for the solver.
  SUBROUTINE FactorizedChebyshevCompute_wrp(ih_InputMat, ih_OutputMat, &
       & ih_polynomial, ih_solver_parameters) &
       & bind(c,name="FactorizedChebyshevCompute_wrp")
    INTEGER(kind=c_int), INTENT(in)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(ChebyshevPolynomial_wrp)     :: h_polynomial
    TYPE(FixedSolverParameters_wrp)   :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_polynomial = TRANSFER(ih_polynomial, h_polynomial)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL FactorizedChebyshevCompute(h_InputMat%data, h_OutputMat%data, &
         & h_polynomial%data, h_solver_parameters%data)
  END SUBROUTINE FactorizedChebyshevCompute_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ChebyshevSolversModule_wrp
