#ifndef SolverParameters_h
#define SolverParameters_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
//!\namespace NTPoly C++ interface to NTPoly
namespace NTPoly {
class Permutation;
class SolverBase;
////////////////////////////////////////////////////////////////////////////////
//! A class to store all the parameters used for solvers.
class SolverParameters {
public:
  //! Constructor.
  SolverParameters();
  //! When do we consider a calculation converged.
  //!\param new_value
  void SetConvergeDiff(double new_value);
  //! Max iterations to perform.
  //!\param new_value
  void SetMaxIterations(int new_value);
  //! Whether to have a verbose calculation.
  //!\param new_value
  void SetVerbosity(bool new_value);
  //! Threshold for flushing small values.
  //!\param new_value
  void SetThreshold(double new_value);
  //! Load balance settings.
  //!\param new_value
  void SetLoadBalance(const Permutation &new_value);
  //! Thresholds for step size searches.
  //!\param new_value
  void SetStepThreshold(double new_value);
  //! Whether to automatically monitor convergence
  //!\param new_value
  void SetMonitorConvergence(bool new_value);
  ~SolverParameters();

private:
  //! Pointer to internal data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  SolverParameters(const SolverParameters &);
  //! Assignment operator, locked.
  SolverParameters &operator=(const SolverParameters &);
  friend class SolverBase;
};
} // namespace NTPoly
#endif
