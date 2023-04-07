#ifndef SQUAREROOT_ch
#define SQUAREROOT_ch

void SquareRoot_wrp(const int *ih_Input, int *ih_Output,
                    const int *ih_solver_parameters);
void DenseSquareRoot_wrp(const int *ih_Input, int *ih_Output,
                         const int *ih_solver_parameters);
void InverseSquareRoot_wrp(const int *ih_Input, int *ih_Output,
                           const int *ih_solver_parameters);
void DenseInverseSquareRoot_wrp(const int *ih_Input, int *ih_Output,
                                const int *ih_solver_parameters);

#endif
