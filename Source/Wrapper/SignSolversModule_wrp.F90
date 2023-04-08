!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the sign solvers module for calling from other languages.
MODULE SignSolversModule_wrp
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SignSolversModule
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SignFunction_wrp
  PUBLIC :: DenseSignFunction_wrp
  PUBLIC :: PolarDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function.
  SUBROUTINE SignFunction_wrp(ih_Mat1, ih_SignMat, ih_solver_parameters) &
       & BIND(c,name="SignFunction_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_SignMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_SignMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_SignMat = TRANSFER(ih_SignMat,h_SignMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SignFunction(h_Mat1%DATA, h_SignMat%DATA, h_solver_parameters%DATA)
  END SUBROUTINE SignFunction_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function (dense versoin).
  SUBROUTINE DenseSignFunction_wrp(ih_Mat1, ih_SignMat, ih_solver_parameters) &
       & BIND(c,name="DenseSignFunction_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_SignMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_SignMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_SignMat = TRANSFER(ih_SignMat,h_SignMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL DenseSignFunction(h_Mat1%DATA, h_SignMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE DenseSignFunction_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the polar decomposition of a matrix Mat1 = U*H.
  SUBROUTINE PolarDecomposition_wrp(ih_Mat1, ih_Umat, ih_Hmat, &
       & ih_solver_parameters) BIND(c,name="PolarDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Umat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Hmat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_Umat
    TYPE(Matrix_ps_wrp) :: h_Hmat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_Umat = TRANSFER(ih_Umat,h_Umat)
    h_Hmat = TRANSFER(ih_Hmat,h_Hmat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PolarDecomposition(h_Mat1%DATA, h_Umat%DATA, h_Hmat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE PolarDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SignSolversModule_wrp
