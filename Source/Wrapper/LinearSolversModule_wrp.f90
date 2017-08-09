!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the matrix inversion module for calling from other languages.
MODULE LinearSolversModule_wrp
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE LinearSolversModule, ONLY : CGSolver
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CGSolver_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Solve the matrix equation AX = B using conjugate gradient.
  !! @param[in] ih_AMat the matrix A, must be symmetric, positive definite.
  !! @param[out] ih_XMat the solved for matrix X.
  !! @param[in] ih_BMat the right hand side.
  !! @param[in] ih_solver_parameters parameters for the solver
  SUBROUTINE CGSolver_wrp(ih_AMat, ih_XMat, ih_BMat, ih_solver_parameters) &
       & bind(c,name="CGSolver_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_AMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_XMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_BMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
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
END MODULE LinearSolversModule_wrp
