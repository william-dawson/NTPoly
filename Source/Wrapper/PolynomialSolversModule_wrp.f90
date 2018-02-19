!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the polynomial solvers module for calling from other languages.
MODULE PolynomialSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE FixedSolversModule_wrp, ONLY : FixedSolverParameters_wrp
  USE PolynomialSolversModule, ONLY : Polynomial_t, ConstructPolynomial, &
       & DestructPolynomial, SetCoefficient, HornerCompute, &
       & PatersonStockmeyerCompute
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
  PURE SUBROUTINE ConstructPolynomial_wrp(ih_this, degree) &
       & bind(c,name="ConstructPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    TYPE(Polynomial_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructPolynomial(h_this%data,degree)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  PURE SUBROUTINE DestructPolynomial_wrp(ih_this) &
       & bind(c,name="DestructPolynomial_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(Polynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructPolynomial(h_this%data)
    DEALLOCATE(h_this%data)
  END SUBROUTINE DestructPolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set coefficient of a polynomial.
  SUBROUTINE SetCoefficient_wrp(ih_this, degree, coefficient) &
       & bind(c,name="SetCoefficient_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    REAL(NTREAL), INTENT(IN) :: coefficient
    TYPE(Polynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetCoefficient(h_this%data, degree, coefficient)
  END SUBROUTINE SetCoefficient_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using Horner's Method.
  SUBROUTINE HornerCompute_wrp(ih_InputMat, ih_OutputMat, ih_polynomial, &
       & ih_solver_parameters) bind(c,name="HornerCompute_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(Polynomial_wrp)              :: h_polynomial
    TYPE(FixedSolverParameters_wrp)   :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_polynomial = TRANSFER(ih_polynomial, h_polynomial)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL HornerCompute(h_InputMat%data, h_OutputMat%data, &
         & h_polynomial%data, h_solver_parameters%data)
  END SUBROUTINE HornerCompute_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Polynomial Using Paterson and Stockmeyer's method.
  SUBROUTINE PatersonStockmeyerCompute_wrp(ih_InputMat, ih_OutputMat, &
       & ih_polynomial, ih_solver_parameters) &
       & bind(c,name="PatersonStockmeyerCompute_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(Polynomial_wrp)              :: h_polynomial
    TYPE(FixedSolverParameters_wrp)   :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_polynomial = TRANSFER(ih_polynomial, h_polynomial)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PatersonStockmeyerCompute(h_InputMat%data, h_OutputMat%data, &
         & h_polynomial%data, h_solver_parameters%data)
  END SUBROUTINE PatersonStockmeyerCompute_wrp
END MODULE PolynomialSolversModule_wrp
