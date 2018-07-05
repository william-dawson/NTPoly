#include "LinearSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "LinearSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void LinearSolvers::CGSolver(
    const Matrix_ps &AMat, Matrix_ps &XMat, const Matrix_ps &BMat,
    const IterativeSolverParameters &solver_parameters) {
  CGSolver_wrp(GetIH(AMat), GetIH(XMat), GetIH(BMat), GetIH(solver_parameters));
}
void LinearSolvers::CholeskyDecomposition(
    const Matrix_ps &AMat, Matrix_ps &LMat,
    const FixedSolverParameters &solver_parameters) {
  CholeskyDecomposition_wrp(GetIH(AMat), GetIH(LMat), GetIH(solver_parameters));
}
void LinearSolvers::PivotedCholeskyDecomposition(
    const Matrix_ps &AMat, Matrix_ps &LMat, int rank,
    const FixedSolverParameters &solver_parameters) {
  PivotedCholeskyDecomposition_wrp(GetIH(AMat), GetIH(LMat), &rank,
                                   GetIH(solver_parameters));
}
} // namespace NTPoly
