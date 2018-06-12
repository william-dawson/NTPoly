!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the Chebyshev Solvers Module
MODULE ChebyshevSolversModule_wrp
  USE ChebyshevSolversModule, ONLY : ChebyshevPolynomial_t, &
       & ConstructChebyshevPolynomial, DestructChebyshevPolynomial, &
       & SetChebyshevCoefficient, ChebyshevCompute, FactorizedChebyshevCompute
  USE MatrixPSModule_wrp, ONLY : Matrix_ps_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE FixedSolversModule_wrp, ONLY : FixedSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the polynomial data type.
  TYPE, PUBLIC :: ChebyshevPolynomial_wrp
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
  PURE SUBROUTINE ConstructChebyshevPolynomial_wrp(ih_this, degree) &
       & bind(c,name="ConstructChebyshevPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    TYPE(ChebyshevPolynomial_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructChebyshevPolynomial(h_this%data,degree)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructChebyshevPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  PURE SUBROUTINE DestructChebyshevPolynomial_wrp(ih_this) &
       & bind(c,name="DestructChebyshevPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(ChebyshevPolynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructChebyshevPolynomial(h_this%data)
    DEALLOCATE(h_this%data)
  END SUBROUTINE DestructChebyshevPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set coefficient of a polynomial.
  SUBROUTINE SetChebyshevCoefficient_wrp(ih_this, degree, coefficient) &
       & bind(c,name="SetChebyshevCoefficient_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    REAL(NTREAL), INTENT(IN) :: coefficient
    TYPE(ChebyshevPolynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetChebyshevCoefficient(h_this%data, degree, coefficient)
  END SUBROUTINE SetChebyshevCoefficient_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Chebyshev Polynomial.
  SUBROUTINE ChebyshevCompute_wrp(ih_InputMat, ih_OutputMat, ih_polynomial, &
       & ih_solver_parameters) bind(c,name="ChebyshevCompute_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
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
  SUBROUTINE FactorizedChebyshevCompute_wrp(ih_InputMat, ih_OutputMat, &
       & ih_polynomial, ih_solver_parameters) &
       & bind(c,name="FactorizedChebyshevCompute_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
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
