#ifndef SIGNSOLVERS_h
#define SIGNSOLVERS_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing The Matrix Sign Function.
class SignSolvers : public SolverBase {
public:
  //! Compute the matrix sign function.
  //!\param Mat1 input matrix.
  //!\param SignMat = Sign(Mat1)
  //!\param solver_parameters parameters for the solver
  static void ComputeSign(const Matrix_ps &Mat1, Matrix_ps &SignMat,
                          const SolverParameters &solver_parameters);
  //! Computes the polar decomposition of a matrix Mat1 = U*H.
  //!\param Mat1 input matrix.
  //!\param Umat the unitary polar factor.
  //!\param Hmat the hermitian matrix factor.
  //!\param solver_parameters parameters for the solver
  static void
  ComputePolarDecomposition(const Matrix_ps &Mat1, Matrix_ps &Umat,
                            Matrix_ps &Hmat,
                            const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
