#ifndef ITERATIVESOLVERPARAMETERS_h
#define ITERATIVESOLVERPARAMETERS_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
//!\namespace NTPoly C++ interface to NTPoly
namespace NTPoly {
class Permutation;
class SolverBase;
////////////////////////////////////////////////////////////////////////////////
//! A class to store all the parameters used for solvers.
class IterativeSolverParameters {
public:
  //! Constructor.
  IterativeSolverParameters();
  //! When do we consider a calculation converged.
  //!\param new_value
  void SetConvergeDiff(double new_value);
  //! Max iterations to perform.
  //!\param new_value
  void SetMaxIterations(int new_value);
  //! Where to have a verbose calculation.
  //!\param new_value
  void SetVerbosity(bool new_value);
  //! Threshold for flushing small values.
  //!\param new_value
  void SetThreshold(double new_value);
  //! Load balance settings.
  //!\param new_value
  void SetLoadBalance(const Permutation &new_value);
  ~IterativeSolverParameters();

private:
  //! Pointer to internal data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  IterativeSolverParameters(const IterativeSolverParameters &);
  //! Assignment operator, locked.
  IterativeSolverParameters &operator=(const IterativeSolverParameters &);
  friend class SolverBase;
};
} // namespace NTPoly
#endif
