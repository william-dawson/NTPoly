#ifndef GEOMETRYOPTIMIZATION_h
#define GEOMETRYOPTIMIZATION_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class Matrix_ps;
//! A Class For Solving Chemistry Systems Based On Sparse Matrices.
class GeometryOptimization : public SolverBase {
public:
  //! Create a new guess at the Density Matrix after updating the geometry.
  //! Based on the purification algorithm in \cite niklasson2010trace .
  //!\param PreviousDensity to extrapolate from.
  //!\param Overlap the overlap matrix of the new geometry.
  //!\param nel the number of electrons.
  //!\param NewDensity the extrapolated density.
  //!\param solver_parameters parameters for the solver
  static void
  PurificationExtrapolate(const Matrix_ps &PreviousDensity,
                          const Matrix_ps &Overlap, int nel,
                          Matrix_ps &NewDensity,
                          const IterativeSolverParameters &solver_parameters);
  //! Create a new guess at the Density Matrix after updating the geometry.
  //! Based on the lowdin algorithm in \cite exner2002comparison .
  //!\param PreviousDensity to extrapolate from.
  //!\param OldOverlap the overlap matrix of the old geometry.
  //!\param NewOverlap the overlap matrix of the new geometry.
  //!\param NewDensity the extrapolated density.
  //!\param solver_parameters parameters for the solver
  static void
  LowdinExtrapolate(const Matrix_ps &PreviousDensity,
                    const Matrix_ps &OldOverlap, const Matrix_ps &NewOverlap,
                    Matrix_ps &NewDensity,
                    const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
