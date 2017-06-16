#ifndef ITERATIVESOLVERPARAMETERS_ch
#define ITERATIVESOLVERPARAMETERS_ch

void ConstructIterativeSolverParameters_wrp(int *ih_this);
void SetIterativeConvergeDiff_wrp(int *ih_this, const double *new_value);
void SetIterativeMaxIterations_wrp(int *ih_this, const int *new_value);
void SetIterativeBeVerbose_wrp(int *ih_this, const bool *new_value);
void SetIterativeThreshold_wrp(int *ih_this, const double *new_value);
void SetIterativeLoadBalance_wrp(int *ih_this, const int *ih_permutation);
void DestructIterativeSolverParameters_wrp(int *ih_this);

#endif
