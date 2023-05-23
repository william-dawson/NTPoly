#ifndef DENSITYMATRIXSOLVERS_ch
#define DENSITYMATRIXSOLVERS_ch

void ComputeDenseFOE_wrp(const int *ih_Hamiltonian,
                         const int *ih_InverseSquareRoot, const double *trace,
                         int *ih_Density, const double *inv_temp_in,
                         const double *energy_value_out,
                         const double *chemical_potential_out,
                         const int *ih_solver_parameters);
void WOM_GC_wrp(const int *ih_Hamiltonian, const int *ih_InverseSquareRoot,
                int *ih_Density, const double *chemical_potential,
                const double *inv_temp, const double *energy_value_out,
                const int *ih_solver_parameters);
void WOM_C_wrp(const int *ih_Hamiltonian, const int *ih_InverseSquareRoot,
               int *ih_Density, const double *trace, const double *inv_temp,
               const double *energy_value_out, const int *ih_solver_parameters);

#endif
