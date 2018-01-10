#include "HermiteSolvers.h"
#include "FixedSolversParameters.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "HermiteSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
HermitePolynomial::HermitePolynomial(int degree) {
  ConstructHermitePolynomial_wrp(this->ih_this, &degree);
}
//////////////////////////////////////////////////////////////////////////////
HermitePolynomial::~HermitePolynomial() {
  DestructHermitePolynomial_wrp(this->ih_this);
}
//////////////////////////////////////////////////////////////////////////////
void HermitePolynomial::SetCoefficient(int degree, double coefficient) {
  int temp_degree = degree + 1;
  SetHermiteCoefficient_wrp(this->ih_this, &temp_degree, &coefficient);
}
//////////////////////////////////////////////////////////////////////////////
void HermitePolynomial::Compute(
    const DistributedSparseMatrix &InputMat, DistributedSparseMatrix &OutputMat,
    const FixedSolverParameters &solver_parameters) const {
  HermiteCompute_wrp(GetIH(InputMat), GetIH(OutputMat), this->ih_this,
                       GetIH(solver_parameters));
}
} // namespace NTPoly
