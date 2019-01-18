#ifndef EIGENBOUNDS_ch
#define EIGENBOUNDS_ch

void GershgorinBounds_wrp(const int *ih_Hamiltonian, double *max_value,
                          double *min_value);
void PowerBounds_wrp(const int *ih_Hamiltonian, double *max_value,
                     const int *ih_solver_parameters);
void InteriorEigenvalues_wrp(const int *ih_mat, const int* ih_density,
                             int* nel, int* nvals, int* ih_vecs,
                             const int *ih_solver_parameters);
void SubspaceIteration_wrp(const int *ih_mat, int* ih_vecs, const int* k,
                           const int *ih_solver_parameters);
#endif
