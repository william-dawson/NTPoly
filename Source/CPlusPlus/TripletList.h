#ifndef TRIPLETLIST_h
#define TRIPLETLIST_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
class DistributedSparseMatrix;
template <class T> class SparseMatrix;
template <class T> class Triplet;
//! A data type for a list of triplets.
//! As this is related to matrix multiplication, the referencing indices are
//! rows and columns.
template <class T> class TripletList {
public:
  //! Construct the triplet list.
  //!\param size the size of the list.
  TripletList(int size = 0);
  //! Construct a triplet list from a distributed sparse matrix.
  //!\param matrix to construct from.
  // TripletList(const DistributedSparseMatrix &matrix);
  //! Increase the size of a triplet list.
  //!\param size the new size.
  void Resize(int size);
  //! Add a value to the end of the triplet list.
  //!\param value the triplet value to append.
  void Append(const Triplet<T> &value);
  //! Set a triplet value.
  //!\param index location to set the triplet at.
  //!\param value the triplet value to set.
  void SetTripletAt(int index, const Triplet<T> &value);
  //! Get the triplet value at a given index.
  //!\param index location to get the triplet at.
  Triplet<T> GetTripletAt(int index) const;
  //! Get the number of entries in a triplet list.
  //!\result the number of entries in the list.
  int GetSize() const;
  //! Standard destructor.
  ~TripletList();
  //! Sort a triplet list
  //!\param list to be sorted.
  //!\param matrix_columns this is the highest column value in the list
  //!\param sorted a now sorted version of the list.
  static void SortTripletList(const TripletList<T> &list, int matrix_columns,
                              TripletList<T> &sorted);

private:
  //! Handle to the actual data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  TripletList(const TripletList<T> &);
  //! Assignment operator, locked.
  TripletList<T> &operator=(const TripletList<T> &);

private:
  template <class T2> friend class SparseMatrix;
  friend class DistributedSparseMatrix;
};
} // namespace NTPoly

#endif
