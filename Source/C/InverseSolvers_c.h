#ifndef INVERSESOLVERS_ch
#define INVERSESOLVERS_ch

void CholeskyInvert_wrp(const int *ih_Matrix, int *ih_Inverse,
                        const int *ih_solver_parameters);
void Invert_wrp(const int *ih_Matrix, int *ih_Inverse,
                const int *ih_solver_parameters);
void PseudoInverse_wrp(const int *ih_Matrix, int *ih_Inverse,
                       const int *ih_solver_parameters);
#endif
