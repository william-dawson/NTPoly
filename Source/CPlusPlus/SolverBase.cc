#include "SolverBase.h"

#include "FixedSolversParameters.h"
#include "IterativeSolversParameters.h"
#include "PSMatrix.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
const int *SolverBase::GetIH(const Matrix_ps &dsm) { return dsm.ih_this; }

////////////////////////////////////////////////////////////////////////////////
int *SolverBase::GetIH(Matrix_ps &dsm) { return dsm.ih_this; }

////////////////////////////////////////////////////////////////////////////////
const int *SolverBase::GetIH(const IterativeSolverParameters &csp) {
  return csp.ih_this;
}

////////////////////////////////////////////////////////////////////////////////
int *SolverBase::GetIH(IterativeSolverParameters &csp) { return csp.ih_this; }

////////////////////////////////////////////////////////////////////////////////
const int *SolverBase::GetIH(const FixedSolverParameters &csp) {
  return csp.ih_this;
}

////////////////////////////////////////////////////////////////////////////////
int *SolverBase::GetIH(FixedSolverParameters &csp) { return csp.ih_this; }
} // namespace NTPoly
