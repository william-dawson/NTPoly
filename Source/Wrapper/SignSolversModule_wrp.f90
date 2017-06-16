!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the overlap solvers module for calling from other languages.
MODULE SignSolversModule_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE SignSolversModule, ONLY : SignFunction
  USE DistributedBlockedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SignFunction_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix.
  !> Computes the matrix sign function.
  !! @param[in] ih_Mat1 the input matrix.
  !! @param[out] ih_SignMat the sign of Mat1.
  !! @param[in]  ih_solver_parameters parameters for the solver.
  SUBROUTINE SignFunction_wrp(ih_Mat1, ih_SignMat, ih_solver_parameters) &
       & bind(c,name="SignFunction_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_SignMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_Mat1
    TYPE(DistributedSparseMatrix_wrp) :: h_SignMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_SignMat = TRANSFER(ih_SignMat,h_SignMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SignFunction(h_Mat1%data, h_SignMat%data, h_solver_parameters%data)
  END SUBROUTINE SignFunction_wrp
END MODULE SignSolversModule_wrp
