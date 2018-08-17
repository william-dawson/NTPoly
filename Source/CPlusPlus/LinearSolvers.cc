#include "LinearSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "LinearSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void LinearSolvers::CGSolver(const Matrix_ps &AMat, Matrix_ps &XMat,
                             const Matrix_ps &BMat,
                             const SolverParameters &solver_parameters) {
  CGSolver_wrp(GetIH(AMat), GetIH(XMat), GetIH(BMat), GetIH(solver_parameters));
}
} // namespace NTPoly
