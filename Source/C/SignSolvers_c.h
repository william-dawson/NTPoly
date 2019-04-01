#ifndef SIGNSOLVERS_ch
#define SIGNSOLVERS_ch

void SignFunction_wrp(const int *ih_mat1, int *ih_signmat,
                      const int *ih_solver_parameters);
void PolarDecomposition_wrp(const int *ih_mat1, int *ih_umat, int *ih_hmat,
                            const int *ih_solver_parameters);
#endif
