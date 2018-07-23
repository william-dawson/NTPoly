#ifndef FERMIOPERATOREXPANSION_ch
#define FERMIOPERATOREXPANSION_ch

void ComputeFOEFixed_wrp(const int *ih_Hamiltonian,
                         const int *ih_InverseSquareRoot, const int *nel,
                         int *ih_Density, const int *degree,
                         const double *chemical_potential,
                         const int *ih_solver_parameters);
void ComputeFOESearch_wrp(const int *ih_Hamiltonian,
                          const int *ih_InverseSquareRoot, const int *nel,
                          int *ih_Density, const int *degree,
                          const int *ih_solver_parameters);

#endif
