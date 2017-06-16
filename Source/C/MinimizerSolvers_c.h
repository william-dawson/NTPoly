#ifndef MINIMIZERSOLVERS_ch
#define MINIMIZERSOLVERS_ch

void ConjugateGradient_wrp(const int *ih_Hamiltonian,
                           const int *ih_InverseSquareRoot, const int *nel,
                           int *ih_Density,
                           const double *chemical_potential_out,
                           const int *ih_solver_parameters);

#endif
