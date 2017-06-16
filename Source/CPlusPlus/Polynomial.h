#ifndef POLYNOMIALSOLVERS_h
#define POLYNOMIALSOLVERS_h

#include "SolverBase.h"
#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class FixedSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing General Matrix Polynomials.
class Polynomial : public SolverBase {
public:
  //! Basic constructor.
  //!\param degree of the polynomial.
  Polynomial(int degree);

public:
  //! Set a polynomial coefficient.
  //!\param degree for which to set the coefficient.
  //!\param coefficient value.
  void SetCoefficient(int degree, double coefficient);

public:
  //! Compute A Matrix Polynomial Using Horner's Method.
  //!\param InputMat input matrix.
  //!\param OutputMat = p(InputMat)
  //!\param solver_parameters parameters for the solver
  void HornerCompute(const DistributedSparseMatrix &InputMat,
                     DistributedSparseMatrix &OutputMat,
                     const FixedSolverParameters &solver_parameters) const;
  //! Compute A Matrix Polynomial Using Paterson and Stockmeyer's Method.
  //!\param InputMat input matrix.
  //!\param OutputMat = p(InputMat)
  //!\param solver_parameters parameters for the solver
  void PatersonStockmeyerCompute(
      const DistributedSparseMatrix &InputMat,
      DistributedSparseMatrix &OutputMat,
      const FixedSolverParameters &solver_parameters) const;

public:
  //! Standard destructor.
  ~Polynomial();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  Polynomial &operator=(const Polynomial &);
  Polynomial(const NTPoly::Polynomial &matB);
};
} // namespace NTPoly
#endif
