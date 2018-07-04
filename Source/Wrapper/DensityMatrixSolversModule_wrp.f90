!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the density matrix solvers module for calling from other languages.
MODULE DensityMatrixSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DensityMatrixSolversModule, ONLY : TRS2, TRS4, HPCP
  USE MatrixPSModule_wrp, ONLY : &
       & Matrix_ps_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TRS2_wrp
  PUBLIC :: TRS4_wrp
  PUBLIC :: HPCP_wrp
  ! PUBLIC :: HPCPPlus_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS2 method.
  SUBROUTINE TRS2_wrp(ih_Hamiltonian, ih_InverseSquareRoot, &
       & nel, ih_Density, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="TRS2_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: chemical_potential_out
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL TRS2(h_Hamiltonian%data, h_InverseSquareRoot%data, INT(nel), &
         & h_Density%data, chemical_potential_out, h_solver_parameters%data)
  END SUBROUTINE TRS2_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS4 method.
  SUBROUTINE TRS4_wrp(ih_Hamiltonian, ih_InverseSquareRoot, &
       & nel, ih_Density, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="TRS4_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: chemical_potential_out
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL TRS4(h_Hamiltonian%data, h_InverseSquareRoot%data, nel, &
         & h_Density%data, chemical_potential_out, h_solver_parameters%data)
  END SUBROUTINE TRS4_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP method.
  SUBROUTINE HPCP_wrp(ih_Hamiltonian, ih_InverseSquareRoot, &
       & nel, ih_Density, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="HPCP_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: chemical_potential_out
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL HPCP(h_Hamiltonian%data, h_InverseSquareRoot%data, nel, &
         & h_Density%data, chemical_potential_out, h_solver_parameters%data)
  END SUBROUTINE HPCP_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DensityMatrixSolversModule_wrp
