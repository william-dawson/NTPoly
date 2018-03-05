!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigensolvers module for calling from other languages.
MODULE EigenSolversModule_wrp
  USE EigenSolversModule, ONLY : EigenDecomposition, SingularValueDecomposition
  USE DistributedSparseMatrixModule_wrp, ONLY : DistributedSparseMatrix_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition_wrp
  PUBLIC :: SingularValueDecompostion_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvalues and eigenvectors of a matrix.
  SUBROUTINE EigenDecomposition_wrp(ih_this, ih_eigenvectors, ih_eigenvalues, &
       & ih_solver_parameters) bind(c,name="EigenDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvalues(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(DistributedSparseMatrix_wrp) :: h_eigenvectors
    TYPE(DistributedSparseMatrix_wrp) :: h_eigenvalues
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_eigenvectors = TRANSFER(ih_eigenvectors,h_eigenvectors)
    h_eigenvalues = TRANSFER(ih_eigenvalues,h_eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL EigenDecomposition(h_this%data, h_eigenvectors%data, &
         & h_eigenvalues%data, h_solver_parameters%data)
  END SUBROUTINE EigenDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the singularvalues and singularvectors of a matrix.
  SUBROUTINE SingularValueDecompostion_wrp(ih_this, ih_leftvectors, &
       & ih_rightvectors, ih_singularvalues, ih_solver_parameters) &
       & bind(c,name="SingularValueDecompostion_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_leftvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_rightvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_singularvalues(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(DistributedSparseMatrix_wrp) :: h_leftvectors
    TYPE(DistributedSparseMatrix_wrp) :: h_rightvectors
    TYPE(DistributedSparseMatrix_wrp) :: h_singularvalues
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_leftvectors = TRANSFER(ih_leftvectors,h_leftvectors)
    h_rightvectors = TRANSFER(ih_rightvectors,h_rightvectors)
    h_singularvalues = TRANSFER(ih_singularvalues,h_singularvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SingularValueDecomposition(h_this%data, h_leftvectors%data, &
         & h_rightvectors%data, h_singularvalues%data, h_solver_parameters%data)
  END SUBROUTINE SingularValueDecompostion_wrp
END MODULE EigenSolversModule_wrp
