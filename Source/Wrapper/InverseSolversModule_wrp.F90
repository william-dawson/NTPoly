!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the matrix inversion module for calling from other languages.
MODULE InverseSolversModule_wrp
  USE InverseSolversModule, ONLY : Invert, DenseInvert, PseudoInverse
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: Invert_wrp
  PUBLIC :: DenseInvert_wrp
  PUBLIC :: PseudoInverse_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix.
  SUBROUTINE Invert_wrp(ih_Mat1, ih_InverseMat, ih_solver_parameters) &
       & BIND(c,name="Invert_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_InverseMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_InverseMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_InverseMat = TRANSFER(ih_InverseMat,h_InverseMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Invert(h_Mat1%DATA, h_InverseMat%DATA, h_solver_parameters%DATA)
  END SUBROUTINE Invert_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix (dense version).
  SUBROUTINE DenseInvert_wrp(ih_Mat1, ih_InverseMat, ih_solver_parameters) &
       & BIND(c,name="DenseInvert_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_InverseMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_InverseMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_InverseMat = TRANSFER(ih_InverseMat,h_InverseMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL DenseInvert(h_Mat1%DATA, h_InverseMat%DATA, h_solver_parameters%DATA)
  END SUBROUTINE DenseInvert_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the pseudoinverse of a matrix.
  SUBROUTINE PseudoInverse_wrp(ih_Mat1, ih_InverseMat, ih_solver_parameters) &
       & BIND(c,name="PseudoInverse_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_InverseMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_InverseMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_InverseMat = TRANSFER(ih_InverseMat,h_InverseMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PseudoInverse(h_Mat1%DATA, h_InverseMat%DATA, h_solver_parameters%DATA)
  END SUBROUTINE PseudoInverse_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE InverseSolversModule_wrp
