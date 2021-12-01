#ifndef INVERSESOLVERS_ch
#define INVERSESOLVERS_ch

void Invert_wrp(const int *ih_Hamiltonian, int *ih_Inverse,
                const int *ih_solver_parameters);
void DenseInvert_wrp(const int *ih_Hamiltonian, int *ih_Inverse,
                     const int *ih_solver_parameters);
void PseudoInverse_wrp(const int *ih_Hamiltonian, int *ih_Inverse,
                       const int *ih_solver_parameters);
#endif
