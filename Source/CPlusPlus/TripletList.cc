#include "Triplet.h"
#include "TripletList.h"
#include <complex.h>

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "TripletList_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
template <> TripletList<double>::TripletList(int size) {
  ConstructTripletList_r_wrp(ih_this, &size);
}
template <> TripletList<double _Complex>::TripletList(int size) {
  ConstructTripletList_c_wrp(ih_this, &size);
}

////////////////////////////////////////////////////////////////////////////////
// TripletList::TripletList(const DistributedSparseMatrix &matrix) {
//   GetMatrixTripletList_wrp(&matrix.ih_this,&ih_this);
// }

////////////////////////////////////////////////////////////////////////////////
template <> void TripletList<double>::Resize(int size) {
  ResizeTripletList_r_wrp(ih_this, &size);
}
template <> void TripletList<double _Complex>::Resize(int size) {
  ResizeTripletList_c_wrp(ih_this, &size);
}

////////////////////////////////////////////////////////////////////////////////
template <> void TripletList<double>::Append(const Triplet<double> &value) {
  AppendToTripletList_r_wrp(ih_this, &(value.index_column), &(value.index_row),
                            &(value.point_value));
}
template <>
void TripletList<double _Complex>::Append(
    const Triplet<double _Complex> &value) {
  AppendToTripletList_c_wrp(ih_this, &(value.index_column), &(value.index_row),
                            &(value.point_value));
}

////////////////////////////////////////////////////////////////////////////////
template <>
void TripletList<double>::SetTripletAt(int index,
                                       const Triplet<double> &value) {
  int adjusted_index = index + 1;
  SetTripletAt_r_wrp(ih_this, &adjusted_index, &(value.index_column),
                     &(value.index_row), &(value.point_value));
}
template <>
void TripletList<double _Complex>::SetTripletAt(
    int index, const Triplet<double _Complex> &value) {
  int adjusted_index = index + 1;
  SetTripletAt_c_wrp(ih_this, &adjusted_index, &(value.index_column),
                     &(value.index_row), &(value.point_value));
}

////////////////////////////////////////////////////////////////////////////////
template <> Triplet<double> TripletList<double>::GetTripletAt(int index) const {
  Triplet<double> temp;
  int adjusted_index = index + 1;
  GetTripletAt_r_wrp(ih_this, &adjusted_index, &(temp.index_column),
                     &(temp.index_row), &(temp.point_value));
  return temp;
}
template <>
Triplet<double _Complex>
TripletList<double _Complex>::GetTripletAt(int index) const {
  Triplet<double _Complex> temp;
  int adjusted_index = index + 1;
  GetTripletAt_c_wrp(ih_this, &adjusted_index, &(temp.index_column),
                     &(temp.index_row), &(temp.point_value));
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
template <> int TripletList<double>::GetSize() const {
  return GetTripletListSize_r_wrp(ih_this);
}
template <> int TripletList<double _Complex>::GetSize() const {
  return GetTripletListSize_c_wrp(ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <> TripletList<double>::~TripletList() {
  DestructTripletList_r_wrp(ih_this);
}
template <> TripletList<double _Complex>::~TripletList() {
  DestructTripletList_c_wrp(ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void TripletList<double>::SortTripletList(const TripletList<double> &input,
                                          int matrix_columns,
                                          TripletList<double> &sorted) {
  SortTripletList_r_wrp(input.ih_this, &matrix_columns, sorted.ih_this);
}
template <>
void TripletList<double _Complex>::SortTripletList(
    const TripletList<double _Complex> &input, int matrix_columns,
    TripletList<double _Complex> &sorted) {
  SortTripletList_c_wrp(input.ih_this, &matrix_columns, sorted.ih_this);
}
} // namespace NTPoly
