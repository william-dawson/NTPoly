#ifndef FIXEDSOLVERPARAMETERS_ch
#define FIXEDSOLVERPARAMETERS_ch

void ConstructFixedSolverParameters_wrp(int *ih_this);
void SetFixedBeVerbose_wrp(int *ih_this, const bool *new_value);
void SetFixedThreshold_wrp(int *ih_this, const double *new_value);
void SetFixedLoadBalance_wrp(int *ih_this, const int *ih_permutation);
void DestructFixedSolverParameters_wrp(int *ih_this);

#endif
