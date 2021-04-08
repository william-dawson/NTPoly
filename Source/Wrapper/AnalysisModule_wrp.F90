!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the linear solvers module for calling from other languages.
MODULE AnalysisModule_wrp
  USE AnalysisModule, ONLY : PivotedCholeskyDecomposition
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PivotedCholeskyDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Symmetric Positive Semi-Definite
  !! matrix.
  SUBROUTINE PivotedCholeskyDecomposition_wrp(ih_AMat, ih_LMat, rank_in, &
       & ih_solver_parameters) bind(c,name="PivotedCholeskyDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_AMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_LMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: rank_in
    TYPE(Matrix_ps_wrp) :: h_AMat
    TYPE(Matrix_ps_wrp) :: h_LMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_AMat = TRANSFER(ih_AMat,h_AMat)
    h_LMat = TRANSFER(ih_LMat,h_LMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PivotedCholeskyDecomposition(h_AMat%data, h_LMat%data, rank_in, &
         & h_solver_parameters%data)
  END SUBROUTINE PivotedCholeskyDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE AnalysisModule_wrp
