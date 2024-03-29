!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the polynomial solvers module for calling from other languages.
MODULE PolynomialSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE PolynomialSolversModule
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the polynomial data type.
  TYPE, PUBLIC :: Polynomial_wrp
     TYPE(Polynomial_t), POINTER :: DATA
  END TYPE Polynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Polynomial type.
  PUBLIC :: ConstructPolynomial_wrp
  PUBLIC :: DestructPolynomial_wrp
  PUBLIC :: SetCoefficient_wrp
  !! Solvers.
  PUBLIC :: HornerCompute_wrp
  PUBLIC :: PatersonStockmeyerCompute_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the empty polynomial constructor.
  SUBROUTINE ConstructPolynomial_wrp(ih_this, degree) &
       & BIND(c,name="ConstructPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    TYPE(Polynomial_wrp) :: h_this

    ALLOCATE(h_this%DATA)
    CALL ConstructPolynomial(h_this%DATA,degree)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  SUBROUTINE DestructPolynomial_wrp(ih_this) &
       & BIND(c,name="DestructPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(Polynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructPolynomial(h_this%DATA)
    DEALLOCATE(h_this%DATA)
  END SUBROUTINE DestructPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set coefficient of a polynomial.
  SUBROUTINE SetCoefficient_wrp(ih_this, degree, coefficient) &
       & BIND(c,name="SetCoefficient_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    REAL(NTREAL), INTENT(IN) :: coefficient
    TYPE(Polynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetCoefficient(h_this%DATA, degree, coefficient)
  END SUBROUTINE SetCoefficient_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using Horner Method.
  SUBROUTINE HornerCompute_wrp(ih_InputMat, ih_OutputMat, ih_polynomial, &
       & ih_solver_parameters) BIND(c,name="HornerCompute_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(Polynomial_wrp)              :: h_polynomial
    TYPE(SolverParameters_wrp)   :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_polynomial = TRANSFER(ih_polynomial, h_polynomial)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Compute(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_polynomial%DATA, h_solver_parameters%DATA)
  END SUBROUTINE HornerCompute_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using Paterson and Stockmeyer method.
  SUBROUTINE PatersonStockmeyerCompute_wrp(ih_InputMat, ih_OutputMat, &
       & ih_polynomial, ih_solver_parameters) &
       & BIND(c,name="PatersonStockmeyerCompute_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(Polynomial_wrp)              :: h_polynomial
    TYPE(SolverParameters_wrp)   :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_polynomial = TRANSFER(ih_polynomial, h_polynomial)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL FactorizedCompute(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_polynomial%DATA, h_solver_parameters%DATA)
  END SUBROUTINE PatersonStockmeyerCompute_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PolynomialSolversModule_wrp
