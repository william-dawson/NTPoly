#ifndef HERMITESOLVERS_ch
#define HERMITESOLVERS_ch

void ConstructHermitePolynomial_wrp(int *ih_polynomial, const int *degree);
void DestructHermitePolynomial_wrp(int *ih_polynomial);
void SetHermiteCoefficient_wrp(int *ih_polynomial, const int *degree,
                                 const double *coefficient);
void HermiteCompute_wrp(const int *ih_InputMat, int *ih_OutputMat,
                          const int *ih_polynomial,
                          const int *ih_solver_parameters);

#endif
