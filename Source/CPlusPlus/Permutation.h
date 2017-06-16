#ifndef PERMUTATION_h
#define PERMUTATION_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class DistributedSparseMatrix;
class IterativeSolverParameters;
class LoadBalancer;

////////////////////////////////////////////////////////////////////////////////
//! A data structure for storing permutations.
class Permutation {
public:
  //! Standard constructor.
  //\param the dimension of the matrix to load balance.
  Permutation(int matrix_dimension);

public:
  //! Fills the load balancer with a default schedule.
  void SetDefaultPermutation();
  //! Fills the load balancer with a reverse schedule.
  void SetReversePermutation();
  //! ills the load balancer with a random schedule.
  void SetRandomPermutation();

public:
  ~Permutation();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];
  //! Stores the size of the matrix to schedule.
  int matrix_dimension;
  //! True if a schedule has been set.
  bool was_filled;

private:
  Permutation(const Permutation &);
  Permutation &operator=(const Permutation &);

private:
  friend class DistributedSparseMatrix;
  friend class IterativeSolverParameters;
  friend class FixedSolverParameters;
  friend class LoadBalancer;
};
} // namespace NTPoly
#endif
