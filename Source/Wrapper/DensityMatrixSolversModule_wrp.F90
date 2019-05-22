!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the density matrix solvers module for calling from other languages.
MODULE DensityMatrixSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DensityMatrixSolversModule, ONLY : TRS2, TRS4, HPCP, PM, ScaleAndFold, &
       & DACPurify, EnergyDensityMatrix
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PM_wrp
  PUBLIC :: TRS2_wrp
  PUBLIC :: TRS4_wrp
  PUBLIC :: HPCP_wrp
  PUBLIC :: ScaleAndFold_wrp
  PUBLIC :: DACPurify_wrp
  PUBLIC :: EnergyDensityMatrix_wrp
  ! PUBLIC :: HPCPPlus_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the PM method.
  SUBROUTINE PM_wrp(ih_Hamiltonian, ih_InverseSquareRoot, nel, ih_Density, &
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="PM_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
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

    CALL PM(h_Hamiltonian%data, h_InverseSquareRoot%data, INT(nel), &
         & h_Density%data, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%data)
  END SUBROUTINE PM_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS2 method.
  SUBROUTINE TRS2_wrp(ih_Hamiltonian, ih_InverseSquareRoot, nel, ih_Density, &
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="TRS2_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
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

    CALL TRS2(h_Hamiltonian%data, h_InverseSquareRoot%data, INT(nel), &
         & h_Density%data, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%data)
  END SUBROUTINE TRS2_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS4 method.
  SUBROUTINE TRS4_wrp(ih_Hamiltonian, ih_InverseSquareRoot, nel, ih_Density, &
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="TRS4_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
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

    CALL TRS4(h_Hamiltonian%data, h_InverseSquareRoot%data, nel, &
         & h_Density%data, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%data)
  END SUBROUTINE TRS4_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP method.
  SUBROUTINE HPCP_wrp(ih_Hamiltonian, ih_InverseSquareRoot, nel, ih_Density, &
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="HPCP_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
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

    CALL HPCP(h_Hamiltonian%data, h_InverseSquareRoot%data, nel, &
         & h_Density%data, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%data)
  END SUBROUTINE HPCP_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the Scale and Fold
  !> method.
  SUBROUTINE ScaleAndFold_wrp(ih_Hamiltonian, ih_InverseSquareRoot, nel, &
       & ih_Density, homo, lumo, energy_value_out, ih_solver_parameters) &
       & bind(c,name="ScaleAndFold_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: energy_value_out
    REAL(NTREAL), INTENT(IN) :: homo
    REAL(NTREAL), INTENT(IN) :: lumo
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ScaleAndFold(h_Hamiltonian%data, h_InverseSquareRoot%data, nel, &
         & h_Density%data, homo, lumo, energy_value_out, &
         & h_solver_parameters%data)
  END SUBROUTINE ScaleAndFold_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix with a divide and conquer approach
  SUBROUTINE DACPurify_wrp(ih_Hamiltonian, ih_InverseSquareRoot, gap_list, &
       & ngaps, ih_Density, energy_value_out, chemical_potential_out, &
       & ih_solver_parameters) bind(c,name="DACPurify_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ngaps
    INTEGER(kind=c_int), DIMENSION(ngaps), INTENT(IN) :: gap_list
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: energy_value_out
    REAL(NTREAL), INTENT(OUT) :: chemical_potential_out
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRoot
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(SolverParameters_wrp) :: h_solver_parameters
    INTEGER(kind=c_int), DIMENSION(ngaps) :: agap_list
    INTEGER :: II

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    DO II = 1, ngaps
       agap_list(II) = gap_list(II)
    END DO

    CALL DACPurify(h_Hamiltonian%data, h_InverseSquareRoot%data, agap_list,&
         & h_Density%data, energy_value_out, chemical_potential_out,&
         & h_solver_parameters%data)
  END SUBROUTINE DACPurify_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the energy-weighted density matrix.
  SUBROUTINE EnergyDensityMatrix_wrp(ih_Hamiltonian, ih_Density, &
       & ih_EnergyDensity, threshold) bind(c,name="EnergyDensityMatrix_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_Density(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_EnergyDensity(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(Matrix_ps_wrp) :: h_Hamiltonian
    TYPE(Matrix_ps_wrp) :: h_Density
    TYPE(Matrix_ps_wrp) :: h_EnergyDensity

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_EnergyDensity = TRANSFER(ih_EnergyDensity,h_EnergyDensity)

    CALL EnergyDensityMatrix(h_Hamiltonian%data, h_Density%data, &
         & h_EnergyDensity%data, threshold)
  END SUBROUTINE EnergyDensityMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DensityMatrixSolversModule_wrp
