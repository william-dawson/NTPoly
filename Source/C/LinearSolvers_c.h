#ifndef LINEARSOLVERS_ch
#define LINEARSOLVERS_ch

void CGSolver_wrp(const int *ih_MatA, int *ih_MatX, const int *ih_matB,
                  const int *ih_solver_parameters);
void CholeskyDecomposition_wrp(const int *ih_MatA, int *ih_MatL,
                               const int *ih_solver_parameters);
#endif
