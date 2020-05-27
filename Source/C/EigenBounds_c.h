#ifndef EIGENBOUNDS_ch
#define EIGENBOUNDS_ch

void GershgorinBounds_wrp(const int *ih_Hamiltonian, double *max_value,
                          double *min_value);
void PowerBounds_wrp(const int *ih_Hamiltonian, double *max_value,
                     const int *ih_solver_parameters);
void DistributedEigenDecomposition_wrp(const int *ih_this, int *ih_eigenvectors,
                                       const int *ih_solver_parameters);
#endif
