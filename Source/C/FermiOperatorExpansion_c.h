#ifndef FERMIOPERATOREXPANSION_ch
#define FERMIOPERATOREXPANSION_ch

void ComputeFOE_wrp(const int *ih_Hamiltonian,
                    const int *ih_InverseSquareRoot, const int *nel,
                    int *ih_Density, const double *energy_value_out,
                    const double *chemical_potential_out,
                    const int *ih_solver_parameters);

#endif
