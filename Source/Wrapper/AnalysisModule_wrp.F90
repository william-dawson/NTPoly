!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the linear solvers module for calling from other languages.
MODULE AnalysisModule_wrp
  USE AnalysisModule, ONLY : PivotedCholeskyDecomposition, ReduceDimension
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PivotedCholeskyDecomposition_wrp
  PUBLIC :: ReduceDimension_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Symmetric Positive Semi-Definite
  !! matrix.
  SUBROUTINE PivotedCholeskyDecomposition_wrp(ih_AMat, ih_LMat, rank_in, &
       & ih_solver_parameters) BIND(c,name="PivotedCholeskyDecomposition_wrp")
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

    CALL PivotedCholeskyDecomposition(h_AMat%DATA, h_LMat%DATA, rank_in, &
         & h_solver_parameters%DATA)
  END SUBROUTINE PivotedCholeskyDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> When we want to only compute the first n eigenvalues of a matrix, this
  !> routine will project out the higher eigenvalues.
  SUBROUTINE ReduceDimension_wrp(ih_this, dim, ih_reduced, &
             & ih_solver_parameters) BIND(c,name="ReduceDimension_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: dim
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_reduced(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_reduced
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_reduced = TRANSFER(ih_reduced,h_reduced)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ReduceDimension(h_this%DATA, dim, h_reduced%DATA, &
          & h_solver_parameters%DATA)
END SUBROUTINE ReduceDimension_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE AnalysisModule_wrp
