#include "TrigonometrySolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "TrigonometrySolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void TrigonometrySolvers::Sine(const Matrix_ps &Input, Matrix_ps &Output,
                               const SolverParameters &solver_parameters) {
  Sine_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void TrigonometrySolvers::Cosine(const Matrix_ps &Input, Matrix_ps &Output,
                                 const SolverParameters &solver_parameters) {
  Cosine_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
} // namespace NTPoly
