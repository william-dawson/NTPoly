#ifndef ChebyshevSOLVERS_ch
#define ChebyshevSOLVERS_ch

void ConstructChebyshevPolynomial_wrp(int *ih_polynomial, const int *degree);
void DestructChebyshevPolynomial_wrp(int *ih_polynomial);
void SetChebyshevCoefficient_wrp(int *ih_polynomial, const int *degree,
                                 const double *coefficient);
void ChebyshevCompute_wrp(const int *ih_InputMat, int *ih_OutputMat,
                          const int *ih_polynomial,
                          const int *ih_solver_parameters);
void FactorizedChebyshevCompute_wrp(const int *ih_InputMat, int *ih_OutputMat,
                                    const int *ih_polynomial,
                                    const int *ih_solver_parameters);

#endif
