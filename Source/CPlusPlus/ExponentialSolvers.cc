#include "ExponentialSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "ExponentialSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void ExponentialSolvers::ComputeExponential(
    const Matrix_ps &Input, Matrix_ps &Output,
    const FixedSolverParameters &solver_parameters) {
  ComputeExponential_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void ExponentialSolvers::ComputeExponentialPade(
    const Matrix_ps &Input, Matrix_ps &Output,
    const IterativeSolverParameters &solver_parameters) {
  ComputeExponentialPade_wrp(GetIH(Input), GetIH(Output),
                             GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void ExponentialSolvers::ComputeLogarithm(
    const Matrix_ps &Input, Matrix_ps &Output,
    const FixedSolverParameters &solver_parameters) {
  ComputeLogarithm_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
} // namespace NTPoly
