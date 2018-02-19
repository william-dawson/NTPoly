!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the sign solvers module for calling from other languages.
MODULE SignSolversModule_wrp
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE SignSolversModule, ONLY : SignFunction, PolarDecomposition
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SignFunction_wrp
  PUBLIC :: PolarDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the matrix sign function.
  SUBROUTINE SignFunction_wrp(ih_Mat1, ih_SignMat, ih_solver_parameters) &
       & bind(c,name="SignFunction_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_SignMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_Mat1
    TYPE(DistributedSparseMatrix_wrp) :: h_SignMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_SignMat = TRANSFER(ih_SignMat,h_SignMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SignFunction(h_Mat1%data, h_SignMat%data, h_solver_parameters%data)
  END SUBROUTINE SignFunction_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes the polar decomposition of a matrix Mat1 = U*H.
  SUBROUTINE PolarDecomposition_wrp(ih_Mat1, ih_Umat, ih_Hmat, &
       & ih_solver_parameters) bind(c,name="PolarDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Umat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Hmat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_Mat1
    TYPE(DistributedSparseMatrix_wrp) :: h_Umat
    TYPE(DistributedSparseMatrix_wrp) :: h_Hmat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_Umat = TRANSFER(ih_Umat,h_Umat)
    h_Hmat = TRANSFER(ih_Hmat,h_Hmat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PolarDecomposition(h_Mat1%data, h_Umat%data, h_Hmat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE PolarDecomposition_wrp
END MODULE SignSolversModule_wrp
