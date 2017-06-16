//! A class to wrap up triplets of values.
#ifndef TRIPLET_h
#define TRIPLET_h

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
//! A Class For Storing Triplets of Integer, Integer, Double.
class Triplet {
public:
  //! Column location.
  int index_column;
  //! Row location.
  int index_row;
  //! Value at that location.
  double point_value;
};

} // namespace NTPoly
#endif
