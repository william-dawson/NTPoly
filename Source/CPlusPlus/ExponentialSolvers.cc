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
    const SolverParameters &solver_parameters) {
  ComputeExponential_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void ExponentialSolvers::ComputeDenseExponential(
    const Matrix_ps &Input, Matrix_ps &Output,
    const SolverParameters &solver_parameters) {
  ComputeDenseExponential_wrp(GetIH(Input), GetIH(Output),
                              GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void ExponentialSolvers::ComputeExponentialPade(
    const Matrix_ps &Input, Matrix_ps &Output,
    const SolverParameters &solver_parameters) {
  ComputeExponentialPade_wrp(GetIH(Input), GetIH(Output),
                             GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void ExponentialSolvers::ComputeLogarithm(
    const Matrix_ps &Input, Matrix_ps &Output,
    const SolverParameters &solver_parameters) {
  ComputeLogarithm_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void ExponentialSolvers::ComputeDenseLogarithm(
    const Matrix_ps &Input, Matrix_ps &Output,
    const SolverParameters &solver_parameters) {
  ComputeDenseLogarithm_wrp(GetIH(Input), GetIH(Output),
                            GetIH(solver_parameters));
}
} // namespace NTPoly
