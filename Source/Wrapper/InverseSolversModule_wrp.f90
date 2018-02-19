!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the matrix inversion module for calling from other languages.
MODULE InverseSolversModule_wrp
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE InverseSolversModule, ONLY : Invert, PseudoInverse
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: Invert_wrp
  PUBLIC :: PseudoInverse_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix.
  SUBROUTINE Invert_wrp(ih_Mat1, ih_InverseMat, ih_solver_parameters) &
       & bind(c,name="Invert_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_InverseMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_Mat1
    TYPE(DistributedSparseMatrix_wrp) :: h_InverseMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_InverseMat = TRANSFER(ih_InverseMat,h_InverseMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Invert(h_Mat1%data, h_InverseMat%data, h_solver_parameters%data)
  END SUBROUTINE Invert_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the pseudoinverse of a matrix.
  SUBROUTINE PseudoInverse_wrp(ih_Mat1, ih_InverseMat, ih_solver_parameters) &
       & bind(c,name="PseudoInverse_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_InverseMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_Mat1
    TYPE(DistributedSparseMatrix_wrp) :: h_InverseMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_InverseMat = TRANSFER(ih_InverseMat,h_InverseMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PseudoInverse(h_Mat1%data, h_InverseMat%data, h_solver_parameters%data)
  END SUBROUTINE PseudoInverse_wrp
END MODULE InverseSolversModule_wrp
