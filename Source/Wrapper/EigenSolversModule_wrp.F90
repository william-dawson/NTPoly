!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigensolvers module for calling from other languages.
MODULE EigenSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE EigenSolversModule
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE PSMAtrixModule, ONLY : PrintMatrix
  USE SingularValueSolversModule, ONLY : SingularValueDecomposition
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : C_INT
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition_wrp
  PUBLIC :: EstimateGap_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigendecomposition of a matrix.
  SUBROUTINE EigenDecomposition_wrp(ih_this, ih_eigenvalues, nvals, &
       & ih_eigenvectors, ih_solver_parameters) &
       & BIND(c,name="EigenDecomposition_wrp")
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_eigenvalues(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: nvals
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_eigenvectors(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_eigenvectors
    TYPE(Matrix_ps_wrp) :: h_eigenvalues
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_eigenvectors = TRANSFER(ih_eigenvectors,h_eigenvectors)
    h_eigenvalues = TRANSFER(ih_eigenvalues,h_eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL EigenDecomposition(h_this%DATA, h_eigenvalues%DATA, nvals_in=nvals, &
         & eigenvectors_in=h_eigenvectors%DATA, &
         & solver_parameters_in=h_solver_parameters%DATA)

  END SUBROUTINE EigenDecomposition_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigendecomposition of a matrix.
  SUBROUTINE EigenDecomposition_novec_wrp(ih_this, ih_eigenvalues, nvals, &
       & ih_solver_parameters) &
       & BIND(c,name="EigenDecomposition_novec_wrp")
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_eigenvalues(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: nvals
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_eigenvalues
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_eigenvalues = TRANSFER(ih_eigenvalues,h_eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL EigenDecomposition(h_this%DATA, h_eigenvalues%DATA, nvals_in=nvals, &
         & solver_parameters_in=h_solver_parameters%DATA)

  END SUBROUTINE EigenDecomposition_novec_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the singularvalues and singularvectors of a matrix.
  SUBROUTINE SingularValueDecompostion_wrp(ih_this, ih_leftvectors, &
       & ih_rightvectors, ih_singularvalues, ih_solver_parameters) &
       & BIND(c,name="SingularValueDecompostion_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_leftvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_rightvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_singularvalues(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_leftvectors
    TYPE(Matrix_ps_wrp) :: h_rightvectors
    TYPE(Matrix_ps_wrp) :: h_singularvalues
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_leftvectors = TRANSFER(ih_leftvectors,h_leftvectors)
    h_rightvectors = TRANSFER(ih_rightvectors,h_rightvectors)
    h_singularvalues = TRANSFER(ih_singularvalues,h_singularvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SingularValueDecomposition(h_this%DATA, h_leftvectors%DATA, &
         & h_rightvectors%DATA, h_singularvalues%DATA, h_solver_parameters%DATA)
  END SUBROUTINE SingularValueDecompostion_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Estimate the HOMO-LUMO gap of a matrix.
  SUBROUTINE EstimateGap_wrp(ih_H, ih_K, chemical_potential, gap, &
       & ih_solver_parameters) &
       & BIND(c,name="EstimateGap_wrp")
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_H(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_K(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    REAL(NTREAL), INTENT(INOUT) :: gap
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_H
    TYPE(Matrix_ps_wrp) :: h_K
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_H = TRANSFER(ih_H, h_H)
    h_K = TRANSFER(ih_K, h_K)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL EstimateGap(h_H%DATA, h_K%DATA, &
         & chemical_potential, gap, h_solver_parameters%DATA)

  END SUBROUTINE EstimateGap_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule_wrp
