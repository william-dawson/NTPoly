#include "RootSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "RootSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
void RootSolvers::ComputeRoot(const Matrix_ps &InputMat, Matrix_ps &OutputMat,
                              int root,
                              const SolverParameters &solver_parameters) {
  ComputeRoot_wrp(GetIH(InputMat), GetIH(OutputMat), &root,
                  GetIH(solver_parameters));
}
//////////////////////////////////////////////////////////////////////////////
void RootSolvers::ComputeInverseRoot(
    const Matrix_ps &InputMat, Matrix_ps &OutputMat, int root,
    const SolverParameters &solver_parameters) {
  ComputeInverseRoot_wrp(GetIH(InputMat), GetIH(OutputMat), &root,
                         GetIH(solver_parameters));
}
} // namespace NTPoly
