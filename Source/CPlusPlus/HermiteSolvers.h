#ifndef HermiteSOLVERS_h
#define HermiteSOLVERS_h

#include "SolverBase.h"
#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing Matrix functions based on Hermite polynomials.
class HermitePolynomial : public SolverBase {
public:
  //! Basic constructor.
  //!\param degree of the polynomial.
  HermitePolynomial(int degree);

public:
  //! Set a polynomial coefficient.
  //!\param degree for which to set the coefficient.
  //!\param coefficient value.
  void SetCoefficient(int degree, double coefficient);

public:
  //! Compute A Matrix Hermite Polynomial.
  //!\param InputMat input matrix.
  //!\param OutputMat = p(InputMat)
  //!\param solver_parameters parameters for the solver
  void Compute(const Matrix_ps &InputMat, Matrix_ps &OutputMat,
               const SolverParameters &solver_parameters) const;

public:
  //! Standard destructor.
  ~HermitePolynomial();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  HermitePolynomial &operator=(const HermitePolynomial &);
  HermitePolynomial(const NTPoly::HermitePolynomial &matB);
};
} // namespace NTPoly
#endif
