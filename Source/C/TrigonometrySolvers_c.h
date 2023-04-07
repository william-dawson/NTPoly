#ifndef TRIGONOMETRY_ch
#define TRIGONOMETRY_ch

void Sine_wrp(const int *ih_Input, int *ih_Output,
              const int *ih_solver_parameters);
void DenseSine_wrp(const int *ih_Input, int *ih_Output,
                   const int *ih_solver_parameters);
void Cosine_wrp(const int *ih_Input, int *ih_Output,
                const int *ih_solver_parameters);
void DenseCosine_wrp(const int *ih_Input, int *ih_Output,
                     const int *ih_solver_parameters);

#endif
