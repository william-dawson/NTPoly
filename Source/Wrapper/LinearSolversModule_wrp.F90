!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the linear solvers module for calling from other languages.
MODULE LinearSolversModule_wrp
  USE LinearSolversModule, ONLY : CGSolver, CholeskyDecomposition
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CGSolver_wrp
  PUBLIC :: CholeskyDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Solve the matrix equation AX = B using conjugate gradient.
  SUBROUTINE CGSolver_wrp(ih_AMat, ih_XMat, ih_BMat, ih_solver_parameters) &
       & BIND(c,name="CGSolver_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_AMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_XMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_BMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_AMat
    TYPE(Matrix_ps_wrp) :: h_XMat
    TYPE(Matrix_ps_wrp) :: h_BMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_AMat = TRANSFER(ih_AMat,h_AMat)
    h_XMat = TRANSFER(ih_XMat,h_XMat)
    h_BMat = TRANSFER(ih_BMat,h_BMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL CGSolver(h_AMat%DATA, h_XMat%DATA, h_BMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE CGSolver_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Symmetric Positive Definite matrix.
  SUBROUTINE CholeskyDecomposition_wrp(ih_AMat, ih_LMat, ih_solver_parameters) &
       & BIND(c,name="CholeskyDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_AMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_LMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_AMat
    TYPE(Matrix_ps_wrp) :: h_LMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_AMat = TRANSFER(ih_AMat,h_AMat)
    h_LMat = TRANSFER(ih_LMat,h_LMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL CholeskyDecomposition(h_AMat%DATA, h_LMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE CholeskyDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LinearSolversModule_wrp
