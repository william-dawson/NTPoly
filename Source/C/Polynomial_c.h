#ifndef POLYNOMIALSOLVERS_ch
#define POLYNOMIALSOLVERS_ch

void ConstructPolynomial_wrp(int *ih_polynomial, const int *degree);
void DestructPolynomial_wrp(int *ih_polynomial);
void SetCoefficient_wrp(int *ih_polynomial, const int *degree,
                        const double *coefficient);
void HornerCompute_wrp(const int *ih_InputMat, int *ih_OutputMat,
                       const int *ih_polynomial,
                       const int *ih_solver_parameters);
void PatersonStockmeyerCompute_wrp(const int *ih_InputMat, int *ih_OutputMat,
                                   const int *ih_polynomial,
                                   const int *ih_solver_parameters);

#endif
