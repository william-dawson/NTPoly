#ifndef ANALYSIS_ch
#define ANALYSIS_ch

void PivotedCholeskyDecomposition_wrp(const int *ih_MatA, int *ih_MatL,
                                      const int *rank_in,
                                      const int *ih_solver_parameters);
void ReduceDimension_wrp(const int *ih_this, const int *dim, int *ih_reduced,
                         const int *ih_solver_parameters);

#endif
