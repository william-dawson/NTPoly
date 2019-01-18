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
void EigenBounds::InteriorEigenvalues(
    const Matrix_ps &matrix, const Matrix_ps &density, int nel, int nvals,
    Matrix_ps &vecs, const SolverParameters &solver_parameters) {
  InteriorEigenvalues_wrp(GetIH(matrix), GetIH(density), &nel, &nvals,
                        GetIH(vecs), GetIH(solver_parameters));
}
void EigenBounds::SubspaceIteration(const Matrix_ps &matrix, Matrix_ps &vecs,
                                    int k,
                                    const SolverParameters &solver_parameters) {
  SubspaceIteration_wrp(GetIH(matrix), GetIH(vecs), &k,
                        GetIH(solver_parameters));
}
} // namespace NTPoly
