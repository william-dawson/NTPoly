//! A class to wrap up triplets of values.
#ifndef TRIPLET_h
#define TRIPLET_h

#include <complex>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
//! A Class For Storing Triplets of Integer, Integer, Double.
class Triplet_r {
public:
  //! Column location.
  int index_column;
  //! Row location.
  int index_row;
  //! Value at that location.
  double point_value;
};

class Triplet_c {
public:
  //! Column location.
  int index_column;
  //! Row location.
  int index_row;
  //! Value at that location.
  std::complex<double> point_value;
};

} // namespace NTPoly
#endif
