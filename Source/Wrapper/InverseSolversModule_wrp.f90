!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the matrix inversion module for calling from other languages.
MODULE InverseSolversModule_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE InverseSolversModule, ONLY : Invert
  USE DistributedBlockedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: Invert_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix.
  !! @param[in]  ih_Mat1 Matrix 1.
  !! @param[out] ih_InverseMat = Mat1^-1.
  !! @param[in]  ih_solver_parameters parameters for the solver
  SUBROUTINE Invert_wrp(ih_Mat1, ih_InverseMat, ih_solver_parameters) &
       & bind(c,name="Invert_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_InverseMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_Mat1
    TYPE(DistributedSparseMatrix_wrp) :: h_InverseMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_InverseMat = TRANSFER(ih_InverseMat,h_InverseMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Invert(h_Mat1%data, h_InverseMat%data, h_solver_parameters%data)
  END SUBROUTINE Invert_wrp
END MODULE InverseSolversModule_wrp
