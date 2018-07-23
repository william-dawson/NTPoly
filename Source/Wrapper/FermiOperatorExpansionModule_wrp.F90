!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the Fermi Operator Expansion Module
MODULE FermiOperatorExpansionModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE FermiOperatorExpansionModule, ONLY : ComputeFOE, FOEEigenvalues
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_double
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeFOEFixed_wrp
  PUBLIC :: ComputeFOESearch_wrp
  PUBLIC :: FOEEigenvalues_wrp
  ! PUBLIC :: FOEEigenvalues_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Density Matrix using the Fermi Operator Expansion.
  SUBROUTINE ComputeFOEFixed_wrp(ih_Hamiltonian, ih_InverseSquareRoot, &
       & nel, ih_Density, degree, chemical_potential, ih_solver_parameters) &
       & bind(c,name="ComputeFOEFixed_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    !! Local
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeFOE(h_Hamiltonian%data, h_InverseSquareRoot%data, INT(nel), &
         & h_Density%data, INT(degree), chemical_potential, &
         & h_solver_parameters%data)
  END SUBROUTINE ComputeFOEFixed_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Density Matrix using the Fermi Operator Expansion.
  SUBROUTINE ComputeFOESearch_wrp(ih_Hamiltonian, ih_InverseSquareRoot, &
       & nel, ih_Density, degree, ih_solver_parameters) &
       & bind(c,name="ComputeFOESearch_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    !! Local
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeFOE(h_Hamiltonian%data, h_InverseSquareRoot%data, INT(nel), &
         & h_Density%data, INT(degree), &
         & solver_parameters_in=h_solver_parameters%data)
  END SUBROUTINE ComputeFOESearch_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Estimate the eigenvalues of a matrix using the Fermi Operator Expansion.
  SUBROUTINE FOEEigenvalues_wrp(ih_InputMat, ih_InverseSquareRoot, &
       & ih_Eigenvalues, degree, nvals, ih_solver_parameters) &
       & bind(c,name="FOEEigenvalues_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_Eigenvalues(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: degree
    INTEGER(kind=c_int), INTENT(IN) :: nvals
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    !! Local
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Eigenvalues
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Eigenvalues = TRANSFER(ih_Eigenvalues,h_Eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL FOEEigenvalues(h_InputMat%data, h_InverseSquareRoot%data, &
         & h_Eigenvalues%data, INT(degree), INT(nvals), &
         & solver_parameters_in=h_solver_parameters%data)
  END SUBROUTINE FOEEigenvalues_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorExpansionModule_wrp
