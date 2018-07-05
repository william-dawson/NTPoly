#ifndef SOLVERBASE_h
#define SOLVERBASE_h

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class Matrix_ps;
class SolverParameters;
////////////////////////////////////////////////////////////////////////////////
//! A class that serves as the base class for all the solvers.
class SolverBase {
protected:
  //! Get the internal handle from a dense sparse matrix.
  //!\param dsm the matrix to get the handle from.
  static const int *GetIH(const Matrix_ps &dsm);
  //! Get the internal handle from a dense sparse matrix.
  //!\param dsm the matrix to get the handle from.
  static int *GetIH(Matrix_ps &dsm);
  //! Get the internal handle from a solver parameter set.
  //!\param csp the parameteres to get the handle from.
  static const int *GetIH(const SolverParameters &csp);
  //! Get the internal handle from a solver parameter set.
  //!\param csp the parameteres to get the handle from.
  static int *GetIH(SolverParameters &csp);
};
} // namespace NTPoly
#endif
