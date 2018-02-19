#ifndef FIXEDSOLVERPARAMETERS_h
#define FIXEDSOLVERPARAMETERS_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
//!\namespace NTPoly C++ interface to NTPoly
namespace NTPoly {
class Permutation;
class SolverBase;
////////////////////////////////////////////////////////////////////////////////
//! A class to store all the parameters used for solvers.
class FixedSolverParameters {
public:
  //! Constructor.
  FixedSolverParameters();
  //! Whether to have a verbose calculation.
  //!\param new_value
  void SetVerbosity(bool new_value);
  //! Threshold for flushing small values.
  //!\param new_value
  void SetThreshold(double new_value);
  //! Load balance settings.
  //!\param new_value
  void SetLoadBalance(const Permutation &new_value);
  ~FixedSolverParameters();

private:
  //! Pointer to internal data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  FixedSolverParameters(const FixedSolverParameters &);
  //! Assignment operator, locked.
  FixedSolverParameters &operator=(const FixedSolverParameters &);
  friend class SolverBase;
};
} // namespace NTPoly
#endif
