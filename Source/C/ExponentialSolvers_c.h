#ifndef EXPONENTIAL_ch
#define EXPONENTIAL_ch

void ComputeExponential_wrp(const int *ih_Input, int *ih_Output,
                            const int *ih_solver_parameters);
void ComputeDenseExponential_wrp(const int *ih_Input, int *ih_Output,
                                 const int *ih_solver_parameters);
void ComputeExponentialPade_wrp(const int *ih_Input, int *ih_Output,
                                const int *ih_solver_parameters);
void ComputeLogarithm_wrp(const int *ih_Input, int *ih_Output,
                          const int *ih_solver_parameters);
void ComputeDenseLogarithm_wrp(const int *ih_Input, int *ih_Output,
                               const int *ih_solver_parameters);

#endif
