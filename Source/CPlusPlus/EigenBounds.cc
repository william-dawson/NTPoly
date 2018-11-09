#include "EigenBounds.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "EigenBounds_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void EigenBounds::GershgorinBounds(const Matrix_ps &matrix, double *min_value,
                                   double *max_value) {
  GershgorinBounds_wrp(GetIH(matrix), min_value, max_value);
}
void EigenBounds::PowerBounds(const Matrix_ps &matrix, double *max_value,
                              const SolverParameters &solver_parameters) {
  PowerBounds_wrp(GetIH(matrix), max_value, GetIH(solver_parameters));
}
} // namespace NTPoly
