!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigensolvers module for calling from other languages.
MODULE EigenSolversModule_wrp
  USE EigenSolversModule, ONLY : EigenDecomposition
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE PSMAtrixModule, ONLY : PrintMatrix
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : C_INT
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigendecomposition of a matrix.
  SUBROUTINE EigenDecomposition_wrp(ih_this, ih_eigenvectors, ih_eigenvalues, &
       & ih_solver_parameters) BIND(c,name="EigenDecomposition_wrp")
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_eigenvectors(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_eigenvalues(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_eigenvectors
    TYPE(Matrix_ps_wrp) :: h_eigenvalues
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_eigenvectors = TRANSFER(ih_eigenvectors,h_eigenvectors)
    h_eigenvalues = TRANSFER(ih_eigenvalues,h_eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL EigenDecomposition(h_this%DATA, h_eigenvectors%DATA, &
         & h_eigenvalues%DATA, h_solver_parameters%DATA)

  END SUBROUTINE EigenDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule_wrp
