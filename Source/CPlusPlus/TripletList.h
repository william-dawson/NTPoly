#ifndef TRIPLETLIST_h
#define TRIPLETLIST_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
class Matrix_ps;
class Triplet_r;
class Triplet_c;
//! A data type for a list of triplets.
//! As this is related to matrix multiplication, the referencing indices are
//! rows and columns.
class TripletList_r {
public:
  //! Construct the triplet list.
  //!\param size the size of the list.
  TripletList_r(int size = 0);
  //! Construct a triplet list from a distributed sparse matrix.
  //!\param matrix to construct from.
  // TripletList(const Matrix_ps &matrix);
  //! Increase the size of a triplet list.
  //!\param size the new size.
  void Resize(int size);
  //! Add a value to the end of the triplet list.
  //!\param value the triplet value to append.
  void Append(const Triplet_r &value);
  //! Set a triplet value.
  //!\param index location to set the triplet at.
  //!\param value the triplet value to set.
  void SetTripletAt(int index, const Triplet_r &value);
  //! Get the triplet value at a given index.
  //!\param index location to get the triplet at.
  Triplet_r GetTripletAt(int index) const;
  //! Get the number of entries in a triplet list.
  //!\result the number of entries in the list.
  int GetSize() const;
  //! Standard destructor.
  ~TripletList_r();
  //! Sort a triplet list
  //!\param list to be sorted.
  //!\param matrix_columns this is the highest column value in the list
  //!\param sorted a now sorted version of the list.
  static void SortTripletList(const TripletList_r &list, int matrix_columns,
                              TripletList_r &sorted);

private:
  //! Handle to the actual data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  TripletList_r(const TripletList_r &);
  //! Assignment operator, locked.
  TripletList_r &operator=(const TripletList_r &);

private:
  friend class Matrix_lsr;
  friend class Matrix_lsc;
  friend class Matrix_ps;
};
//! A data type for a list of triplets.
//! As this is related to matrix multiplication, the referencing indices are
//! rows and columns.
class TripletList_c {
public:
  //! Construct the triplet list.
  //!\param size the size of the list.
  TripletList_c(int size = 0);
  //! Construct a triplet list from a distributed sparse matrix.
  //!\param matrix to construct from.
  // TripletList(const Matrix_ps &matrix);
  //! Increase the size of a triplet list.
  //!\param size the new size.
  void Resize(int size);
  //! Add a value to the end of the triplet list.
  //!\param value the triplet value to append.
  void Append(const Triplet_c &value);
  //! Set a triplet value.
  //!\param index location to set the triplet at.
  //!\param value the triplet value to set.
  void SetTripletAt(int index, const Triplet_c &value);
  //! Get the triplet value at a given index.
  //!\param index location to get the triplet at.
  Triplet_c GetTripletAt(int index) const;
  //! Get the number of entries in a triplet list.
  //!\result the number of entries in the list.
  int GetSize() const;
  //! Standard destructor.
  ~TripletList_c();
  //! Sort a triplet list
  //!\param list to be sorted.
  //!\param matrix_columns this is the highest column value in the list
  //!\param sorted a now sorted version of the list.
  static void SortTripletList(const TripletList_c &list, int matrix_columns,
                              TripletList_c &sorted);

private:
  //! Handle to the actual data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  TripletList_c(const TripletList_c &);
  //! Assignment operator, locked.
  TripletList_c &operator=(const TripletList_c &);

private:
  friend class Matrix_lsr;
  friend class Matrix_lsc;
  friend class Matrix_ps;
};
} // namespace NTPoly

#endif
