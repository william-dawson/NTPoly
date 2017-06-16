#ifndef ROOTSOLVERS_ch
#define ROOTSOLVERS_ch

void ComputeRoot_wrp(const int *ih_inputmat, int *ih_outputmat, const int *root,
                     const int *ih_solver_parameters);
void ComputeInverseRoot_wrp(const int *ih_inputmat, int *ih_outputmat,
                            const int *root, const int *ih_solver_parameters);

#endif
