#ifndef DENSITYMATRIXSOLVERS_ch
#define DENSITYMATRIXSOLVERS_ch

void PM_wrp(const int *ih_Hamiltonian, const int *ih_InverseSquareRoot,
            const double *trace, int *ih_Density,
            const double *energy_value_out,
            const double *chemical_potential_out,
            const int *ih_solver_parameters);
void TRS2_wrp(const int *ih_Hamiltonian, const int *ih_InverseSquareRoot,
              const double *trace, int *ih_Density,
              const double *energy_value_out,
              const double *chemical_potential_out,
              const int *ih_solver_parameters);
void TRS4_wrp(const int *ih_Hamiltonian, const int *ih_InverseSquareRoot,
              const double *trace, int *ih_Density,
              const double *energy_value_out,
              const double *chemical_potential_out,
              const int *ih_solver_parameters);
void HPCP_wrp(const int *ih_Hamiltonian, const int *ih_InverseSquareRoot,
              const double *trace, int *ih_Density,
              const double *energy_value_out,
              const double *chemical_potential_out,
              const int *ih_solver_parameters);
void ScaleAndFold_wrp(const int *ih_Hamiltonian,
                      const int *ih_InverseSquareRoot, const double *trace,
                      int *ih_Density, const double *homo, const double *lumo,
                      const double *energy_value_out,
                      const int *ih_solver_parameters);
void DenseDensity_wrp(const int *ih_Hamiltonian,
                      const int *ih_InverseSquareRoot, const double *trace,
                      int *ih_Density, const double *energy_value_out,
                      const double *chemical_potential_out,
                      const int *ih_solver_parameters);
void EnergyDensityMatrix_wrp(const int *ih_Hamiltonian, const int *ih_Density,
                             int *ih_EnergyDensity, const double *threshold);
void McWeenyStep_wrp(const int *ih_D, int *ih_DOut, const double *threshold);
void McWeenyStepS_wrp(const int *ih_D, int *ih_DOut, const int *ih_S,
                      const double *threshold);
#endif
