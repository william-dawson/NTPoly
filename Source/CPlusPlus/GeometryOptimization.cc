#include "GeometryOptimization.h"
#include "DistributedSparseMatrix.h"
#include "IterativeSolversParameters.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "GeometryOptimization_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void GeometryOptimization::ExtrapolateGeometry(
    const DistributedSparseMatrix &PreviousDensity,
    const DistributedSparseMatrix &Overlap, int nel,
    DistributedSparseMatrix &NewDensity,
    const IterativeSolverParameters &solver_parameters) {
  ExtrapolateGeometry_wrp(GetIH(PreviousDensity), GetIH(Overlap), &nel,
                          GetIH(NewDensity), GetIH(solver_parameters));
}
} // namespace NTPoly
