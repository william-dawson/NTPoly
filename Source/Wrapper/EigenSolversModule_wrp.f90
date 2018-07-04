!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigen solvers module for calling from other languages.
MODULE EigenSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixPSModule_wrp, ONLY : &
       & Matrix_ps_wrp
  USE EigenSolversModule, ONLY : ReferenceEigenDecomposition, &
       & SplittingEigenDecomposition, SingularValueDecomposition
  USE FixedSolversModule_wrp, ONLY : FixedSolverParameters_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SplittingEigenDecomposition_wrp
  PUBLIC :: ReferenceEigenDecomposition_wrp
  PUBLIC :: SingularValueDecompostion_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a distributed matrix.
  SUBROUTINE SplittingEigenDecomposition_wrp(ih_this, ih_eigenvectors, &
       & ih_eigenvalues, target_values, ih_solver_parameters) &
       & bind(c,name="SplittingEigenDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvalues(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: target_values
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_eigenvectors
    TYPE(Matrix_ps_wrp) :: h_eigenvalues
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_eigenvectors = TRANSFER(ih_eigenvectors,h_eigenvectors)
    h_eigenvalues = TRANSFER(ih_eigenvalues,h_eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SplittingEigenDecomposition(h_this%data, h_eigenvectors%data, &
         & h_eigenvalues%data, target_values, h_solver_parameters%data)
  END SUBROUTINE SplittingEigenDecomposition_wrp
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
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_leftvectors
    TYPE(Matrix_ps_wrp) :: h_rightvectors
    TYPE(Matrix_ps_wrp) :: h_singularvalues
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_leftvectors = TRANSFER(ih_leftvectors,h_leftvectors)
    h_rightvectors = TRANSFER(ih_rightvectors,h_rightvectors)
    h_singularvalues = TRANSFER(ih_singularvalues,h_singularvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SingularValueDecomposition(h_this%data, h_leftvectors%data, &
         & h_rightvectors%data, h_singularvalues%data, h_solver_parameters%data)
  END SUBROUTINE SingularValueDecompostion_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a distributed matrix.
  SUBROUTINE ReferenceEigenDecomposition_wrp(ih_this, ih_eigenvectors, &
       & ih_eigenvalues, ih_solver_parameters) &
       & bind(c,name="ReferenceEigenDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvalues(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_eigenvectors
    TYPE(Matrix_ps_wrp) :: h_eigenvalues
    TYPE(FixedSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this, h_this)
    h_eigenvectors = TRANSFER(ih_eigenvectors, h_eigenvectors)
    h_eigenvalues = TRANSFER(ih_eigenvalues, h_eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ReferenceEigenDecomposition(h_this%data, h_eigenvectors%data, &
         & h_eigenvalues%data, h_solver_parameters%data)
  END SUBROUTINE ReferenceEigenDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule_wrp
