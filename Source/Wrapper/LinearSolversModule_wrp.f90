!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the linear solvers module for calling from other languages.
MODULE LinearSolversModule_wrp
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE FixedSolversModule_wrp, ONLY : FixedSolverParameters_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE LinearSolversModule, ONLY : CGSolver, CholeskyDecomposition, &
       & PivotedCholeskyDecomposition
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CGSolver_wrp
  PUBLIC :: CholeskyDecomposition_wrp
  PUBLIC :: PivotedCholeskyDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Solve the matrix equation AX = B using conjugate gradient.
  SUBROUTINE CGSolver_wrp(ih_AMat, ih_XMat, ih_BMat, ih_solver_parameters) &
       & bind(c,name="CGSolver_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_AMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_XMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_BMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_AMat
    TYPE(DistributedSparseMatrix_wrp) :: h_XMat
    TYPE(DistributedSparseMatrix_wrp) :: h_BMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_AMat = TRANSFER(ih_AMat,h_AMat)
    h_XMat = TRANSFER(ih_XMat,h_XMat)
    h_BMat = TRANSFER(ih_BMat,h_BMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL CGSolver(h_AMat%data, h_XMat%data, h_BMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE CGSolver_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Symmetric Positive Definite matrix.
  SUBROUTINE CholeskyDecomposition_wrp(ih_AMat, ih_LMat, ih_solver_parameters) &
       & bind(c,name="CholeskyDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_AMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_LMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_AMat
    TYPE(DistributedSparseMatrix_wrp) :: h_LMat
    TYPE(FixedSolverParameters_wrp) :: h_solver_parameters

    h_AMat = TRANSFER(ih_AMat,h_AMat)
    h_LMat = TRANSFER(ih_LMat,h_LMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL CholeskyDecomposition(h_AMat%data, h_LMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE CholeskyDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Cholesky Decomposition of a Symmetric Positive Semi-Definite matrix.
  SUBROUTINE PivotedCholeskyDecomposition_wrp(ih_AMat, ih_LMat, rank_in, &
       & ih_solver_parameters) bind(c,name="PivotedCholeskyDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_AMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_LMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: rank_in
    TYPE(DistributedSparseMatrix_wrp) :: h_AMat
    TYPE(DistributedSparseMatrix_wrp) :: h_LMat
    TYPE(FixedSolverParameters_wrp) :: h_solver_parameters

    h_AMat = TRANSFER(ih_AMat,h_AMat)
    h_LMat = TRANSFER(ih_LMat,h_LMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PivotedCholeskyDecomposition(h_AMat%data, h_LMat%data, rank_in, &
         & h_solver_parameters%data)
  END SUBROUTINE PivotedCholeskyDecomposition_wrp
END MODULE LinearSolversModule_wrp
