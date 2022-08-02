#include "Analysis.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "Analysis_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
void Analysis::PivotedCholeskyDecomposition(
    const Matrix_ps &AMat, Matrix_ps &LMat, int rank,
    const SolverParameters &solver_parameters) {
  PivotedCholeskyDecomposition_wrp(GetIH(AMat), GetIH(LMat), &rank,
                                   GetIH(solver_parameters));
}
void Analysis::ReduceDimension(const Matrix_ps &AMat, int dim, Matrix_ps &RMat,
    const SolverParameters &solver_parameters) {
  ReduceDimension_wrp(GetIH(AMat), &dim, GetIH(RMat),
                            GetIH(solver_parameters));
}
} // namespace NTPoly
