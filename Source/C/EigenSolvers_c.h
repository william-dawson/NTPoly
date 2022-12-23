#ifndef EIGENSOLVERS_ch
#define EIGENSOLVERS_ch

void EigenDecomposition_wrp(const int *ih_this, int *ih_eigenvectors,
                            const int *nvals, int *ih_eigenvalues,
                            const int *ih_solver_parameters);
void EigenDecomposition_novec_wrp(const int *ih_this, int *ih_eigenvectors,
                                  const int *nvals,
                                  const int *ih_solver_parameters);
void SingularValueDecompostion_wrp(const int *ih_this, int *ih_leftvectors,
                                   int *ih_rightvectors, int *ih_singularvalues,
                                   const int *ih_solver_parameters);
#endif
