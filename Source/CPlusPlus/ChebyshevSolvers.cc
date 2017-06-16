#include "ChebyshevSolvers.h"
#include "FixedSolversParameters.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "ChebyshevSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
ChebyshevPolynomial::ChebyshevPolynomial(int degree) {
  ConstructChebyshevPolynomial_wrp(this->ih_this, &degree);
}
//////////////////////////////////////////////////////////////////////////////
ChebyshevPolynomial::~ChebyshevPolynomial() {
  DestructChebyshevPolynomial_wrp(this->ih_this);
}
//////////////////////////////////////////////////////////////////////////////
void ChebyshevPolynomial::SetCoefficient(int degree, double coefficient) {
  int temp_degree = degree + 1;
  SetChebyshevCoefficient_wrp(this->ih_this, &temp_degree, &coefficient);
}
//////////////////////////////////////////////////////////////////////////////
void ChebyshevPolynomial::Compute(
    const DistributedSparseMatrix &InputMat, DistributedSparseMatrix &OutputMat,
    const FixedSolverParameters &solver_parameters) const {
  ChebyshevCompute_wrp(GetIH(InputMat), GetIH(OutputMat), this->ih_this,
                       GetIH(solver_parameters));
}
//////////////////////////////////////////////////////////////////////////////
void ChebyshevPolynomial::ComputeFactorized(
    const DistributedSparseMatrix &InputMat, DistributedSparseMatrix &OutputMat,
    const FixedSolverParameters &solver_parameters) const {
  FactorizedChebyshevCompute_wrp(GetIH(InputMat), GetIH(OutputMat),
                                 this->ih_this, GetIH(solver_parameters));
}
} // namespace NTPoly
