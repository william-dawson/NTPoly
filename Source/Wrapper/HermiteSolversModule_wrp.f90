!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the Hermite Solvers Module
MODULE HermiteSolversModule_wrp
  USE HermiteSolversModule, ONLY : HermitePolynomial_t, &
       & ConstructHermitePolynomial, DestructHermitePolynomial, &
       & SetHermiteCoefficient, HermiteCompute
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
  TYPE, PUBLIC :: HermitePolynomial_wrp
     !> Actual data.
     TYPE(HermitePolynomial_t), POINTER :: DATA
  END TYPE HermitePolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Hermite Polynomial type.
  PUBLIC :: ConstructHermitePolynomial_wrp
  PUBLIC :: DestructHermitePolynomial_wrp
  PUBLIC :: SetHermiteCoefficient_wrp
  !! Solvers.
  PUBLIC :: HermiteCompute_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the empty polynomial constructor.
  !! @param[out] ih_this handle to the polynomial being created.
  !! @param[in] degree the degree of the polynomial.
  PURE SUBROUTINE ConstructHermitePolynomial_wrp(ih_this, degree) &
       & bind(c,name="ConstructHermitePolynomial_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: degree
    TYPE(HermitePolynomial_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructHermitePolynomial(h_this%data,degree)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructHermitePolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a polynomial object.
  !! @param[inout] ih_this handle to the polynomial to free up.
  PURE SUBROUTINE DestructHermitePolynomial_wrp(ih_this) &
       & bind(c,name="DestructHermitePolynomial_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(HermitePolynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructHermitePolynomial(h_this%data)
    DEALLOCATE(h_this%data)
  END SUBROUTINE DestructHermitePolynomial_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set coefficient of a polynomial.
  !! @param[inout] ih_this handle to the polynomial to set.
  !! @param[in] degree for which to set the coefficient.
  !! @param[in] coefficient value.
  SUBROUTINE SetHermiteCoefficient_wrp(ih_this, degree, coefficient) &
       & bind(c,name="SetHermiteCoefficient_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: degree
    REAL(NTREAL), INTENT(in) :: coefficient
    TYPE(HermitePolynomial_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetHermiteCoefficient(h_this%data, degree, coefficient)
  END SUBROUTINE SetHermiteCoefficient_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute A Matrix Hermite Polynomial.
  !! @param[in] ih_InputMat the input matrix
  !! @param[out] ih_OutputMat = poly(InputMat)
  !! @param[in] ih_polynomial polynomial to compute.
  !! @param[in] ih_solver_parameters parameters for the solver.
  SUBROUTINE HermiteCompute_wrp(ih_InputMat, ih_OutputMat, ih_polynomial, &
       & ih_solver_parameters) bind(c,name="HermiteCompute_wrp")
    INTEGER(kind=c_int), INTENT(in)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_polynomial(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(HermitePolynomial_wrp)     :: h_polynomial
    TYPE(FixedSolverParameters_wrp)   :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_polynomial = TRANSFER(ih_polynomial, h_polynomial)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL HermiteCompute(h_InputMat%data, h_OutputMat%data, &
         & h_polynomial%data, h_solver_parameters%data)
  END SUBROUTINE HermiteCompute_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE HermiteSolversModule_wrp
