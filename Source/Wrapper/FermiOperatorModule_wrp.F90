!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the density matrix solvers module for calling from other languages.
MODULE FermiOperatorModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE FermiOperatorModule
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeDenseFOE_wrp
  PUBLIC :: WOM_GC_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the PM method.
  SUBROUTINE ComputeDenseFOE_wrp(ih_Hamiltonian, ih_InverseSquareRoot, trace, &
       & ih_Density, inv_temp_in, energy_value_out, chemical_potential_out, &
       & ih_solver_parameters) BIND(c,name="ComputeDenseFOE_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: trace
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: inv_temp_in
    REAL(NTREAL), INTENT(OUT) :: energy_value_out
    REAL(NTREAL), INTENT(OUT) :: chemical_potential_out
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeDenseFOE(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, &
         & trace, h_Density%DATA, inv_temp_in=inv_temp_in, &
         & energy_value_out=energy_value_out, &
         & chemical_potential_out=chemical_potential_out, &
         & solver_parameters_in=h_solver_parameters%DATA)
  END SUBROUTINE ComputeDenseFOE_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE WOM_GC_wrp(ih_Hamiltonian, ih_InverseSquareRoot, &
       & ih_Density, chemical_potential, inv_temp, energy_value_out, &
       & ih_solver_parameters) BIND(c,name="WOM_GC_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    REAL(NTREAL), INTENT(IN) :: inv_temp
    REAL(NTREAL), INTENT(OUT) :: energy_value_out
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL WOM_GC(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, &
         & h_Density%DATA, chemical_potential, inv_temp, &
         & energy_value_out=energy_value_out, &
         & solver_parameters_in=h_solver_parameters%DATA)
  END SUBROUTINE WOM_GC_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorModule_wrp
