#ifndef ITERATIVESOLVERPARAMETERS_ch
#define ITERATIVESOLVERPARAMETERS_ch

void ConstructSolverParameters_wrp(int *ih_this);
void SetParametersConvergeDiff_wrp(int *ih_this, const double *new_value);
void SetParametersMaxIterations_wrp(int *ih_this, const int *new_value);
void SetParametersBeVerbose_wrp(int *ih_this, const bool *new_value);
void SetParametersThreshold_wrp(int *ih_this, const double *new_value);
void SetParametersDACBaseSize_wrp(int *ih_this, const int *new_value);
void SetParametersDACBaseSparsity_wrp(int *ih_this, const double *new_value);
void SetParametersLoadBalance_wrp(int *ih_this, const int *ih_permutation);
void DestructSolverParameters_wrp(int *ih_this);

#endif
