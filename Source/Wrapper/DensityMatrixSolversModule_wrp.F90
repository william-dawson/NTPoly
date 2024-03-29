!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the density matrix solvers module for calling from other languages.
MODULE DensityMatrixSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DensityMatrixSolversModule
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
  PUBLIC :: DenseDensity_wrp
  PUBLIC :: ScaleAndFold_wrp
  PUBLIC :: EnergyDensityMatrix_wrp
  PUBLIC :: McWeenyStep_wrp
  ! PUBLIC :: HPCPPlus_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the PM method.
  SUBROUTINE PM_wrp(ih_Hamiltonian, ih_InverseSquareRoot, trace, ih_Density, &
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & BIND(c,name="PM_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: trace
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

    CALL PM(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, trace, &
         & h_Density%DATA, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%DATA)
  END SUBROUTINE PM_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS2 method.
  SUBROUTINE TRS2_wrp(ih_Hamiltonian, ih_InverseSquareRoot, trace, ih_Density,&
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & BIND(c,name="TRS2_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: trace
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

    CALL TRS2(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, trace, &
         & h_Density%DATA, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%DATA)
  END SUBROUTINE TRS2_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS4 method.
  SUBROUTINE TRS4_wrp(ih_Hamiltonian, ih_InverseSquareRoot, trace, ih_Density, &
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & BIND(c,name="TRS4_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: trace
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

    CALL TRS4(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, trace, &
         & h_Density%DATA, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%DATA)
  END SUBROUTINE TRS4_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP method.
  SUBROUTINE HPCP_wrp(ih_Hamiltonian, ih_InverseSquareRoot, trace, ih_Density, &
       & energy_value_out, chemical_potential_out, ih_solver_parameters) &
       & BIND(c,name="HPCP_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: trace
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

    CALL HPCP(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, trace, &
         & h_Density%DATA, energy_value_out, chemical_potential_out, &
         & h_solver_parameters%DATA)
  END SUBROUTINE HPCP_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the Scale and Fold
  !> method.
  SUBROUTINE ScaleAndFold_wrp(ih_Hamiltonian, ih_InverseSquareRoot, trace, &
       & ih_Density, homo, lumo, energy_value_out, ih_solver_parameters) &
       & BIND(c,name="ScaleAndFold_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: trace
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

    CALL ScaleAndFold(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, trace, &
         & h_Density%DATA, homo, lumo, energy_value_out, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ScaleAndFold_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the PM method.
  SUBROUTINE DenseDensity_wrp(ih_Hamiltonian, ih_InverseSquareRoot, trace, &
       & ih_Density, energy_value_out, chemical_potential_out, &
       & ih_solver_parameters) BIND(c,name="DenseDensity_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: trace
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

    CALL DenseDensity(h_Hamiltonian%DATA, h_InverseSquareRoot%DATA, trace, &
         & h_Density%DATA, energy_value_out=energy_value_out, &
         & chemical_potential_out=chemical_potential_out, &
         & solver_parameters_in=h_solver_parameters%DATA)
  END SUBROUTINE DenseDensity_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the energy-weighted density matrix.
  SUBROUTINE EnergyDensityMatrix_wrp(ih_Hamiltonian, ih_Density, &
       & ih_EnergyDensity, threshold) BIND(c,name="EnergyDensityMatrix_wrp")
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

    CALL EnergyDensityMatrix(h_Hamiltonian%DATA, h_Density%DATA, &
         & h_EnergyDensity%DATA, threshold)
  END SUBROUTINE EnergyDensityMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take one McWeeny Step DOut = 3*DD - 2*DDD
  SUBROUTINE McWeenyStep_wrp(ih_D, ih_DOut, threshold) &
       & BIND(c,name="McWeenyStep_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_D(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_DOut(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(Matrix_ps_wrp) :: h_D
    TYPE(Matrix_ps_wrp) :: h_DOut

    h_D = TRANSFER(ih_D,h_D)
    h_DOut = TRANSFER(ih_DOut,h_DOut)

    CALL McWeenyStep(h_D%DATA, h_Dout%DATA, threshold_in=threshold)
  END SUBROUTINE McWeenyStep_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take one McWeeny Step DOut = 3*DSD - 2*DSDSD
  SUBROUTINE McWeenyStepS_wrp(ih_D, ih_DOut, ih_S, threshold) &
       & BIND(c,name="McWeenyStepS_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_D(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_DOut(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_S(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(Matrix_ps_wrp) :: h_D, h_S
    TYPE(Matrix_ps_wrp) :: h_DOut

    h_D = TRANSFER(ih_D,h_D)
    h_S = TRANSFER(ih_S,h_S)
    h_DOut = TRANSFER(ih_DOut,h_DOut)

    CALL McWeenyStep(h_D%DATA, h_Dout%DATA, S_in=h_S%DATA, &
         & threshold_in=threshold)
  END SUBROUTINE McWeenyStepS_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DensityMatrixSolversModule_wrp
