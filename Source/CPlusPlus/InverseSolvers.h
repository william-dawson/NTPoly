#ifndef INVERSESOLVERS_h
#define INVERSESOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A class for inverting matrices.
class InverseSolvers : public SolverBase {
public:
  //! Compute the inverse of a matrix.
  //! An implementation of Hotelling method.
  //!\param Overlap the matrix to invert.
  //!\param InverseMat = Overlap^-1.
  //!\param solver_parameters parameters for the solver
  static void Invert(const Matrix_ps &Overlap, Matrix_ps &InverseMat,
                     const SolverParameters &solver_parameters);
  //! Compute the inverse of a matrix (dense version).
  //!\param Overlap the matrix to invert.
  //!\param InverseMat = Overlap^-1.
  //!\param solver_parameters parameters for the solver
  static void DenseInvert(const Matrix_ps &Overlap, Matrix_ps &InverseMat,
                          const SolverParameters &solver_parameters);
  //! Compute the pseudoinverse of a matrix.
  //! An implementation of Hotelling method, with a different convergence
  //! criteria.
  //!\param Overlap the matrix to invert.
  //!\param InverseMat = Overlap^-1.
  //!\param solver_parameters parameters for the solver
  static void PseudoInverse(const Matrix_ps &Overlap, Matrix_ps &InverseMat,
                            const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
