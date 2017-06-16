#ifndef ChebyshevSOLVERS_h
#define ChebyshevSOLVERS_h

#include "SolverBase.h"
#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class FixedSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing Matrix functions based on Chebyshev polynomials.
class ChebyshevPolynomial : public SolverBase {
public:
  //! Basic constructor.
  //!\param degree of the polynomial.
  ChebyshevPolynomial(int degree);

public:
  //! Set a polynomial coefficient.
  //!\param degree for which to set the coefficient.
  //!\param coefficient value.
  void SetCoefficient(int degree, double coefficient);

public:
  //! Compute A Matrix Chebyshev Polynomial.
  //!\param InputMat input matrix.
  //!\param OutputMat = p(InputMat)
  //!\param solver_parameters parameters for the solver
  void Compute(const DistributedSparseMatrix &InputMat,
               DistributedSparseMatrix &OutputMat,
               const FixedSolverParameters &solver_parameters) const;
  //! Compute A Matrix Chebyshev Polynomial Recursively.
  //!\param InputMat input matrix.
  //!\param OutputMat = p(InputMat)
  //!\param solver_parameters parameters for the solver
  void ComputeFactorized(const DistributedSparseMatrix &InputMat,
                         DistributedSparseMatrix &OutputMat,
                         const FixedSolverParameters &solver_parameters) const;

public:
  //! Standard destructor.
  ~ChebyshevPolynomial();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  ChebyshevPolynomial &operator=(const ChebyshevPolynomial &);
  ChebyshevPolynomial(const NTPoly::ChebyshevPolynomial &matB);
};
} // namespace NTPoly
#endif
