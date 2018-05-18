!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigen solvers module for calling from other languages.
MODULE EigenSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE EigenSolversModule, ONLY : DistributedEigenDecomposition
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: DistributedEigenDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a distributed matrix.
  SUBROUTINE DistributedEigenDecomposition_wrp(ih_this, ih_eigenvectors, &
       & ih_solver_parameters) bind(c,name="DistributedEigenDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(DistributedSparseMatrix_wrp) :: h_eigenvectors
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this, h_this)
    h_eigenvectors = TRANSFER(ih_eigenvectors, h_eigenvectors)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL DistributedEigenDecomposition(h_this%data, h_eigenvectors%data, &
         & h_solver_parameters%data)
  END SUBROUTINE DistributedEigenDecomposition_wrp
END MODULE EigenSolversModule_wrp
