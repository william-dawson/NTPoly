#include "Triplet.h"
#include "TripletList.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "TripletList_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
TripletList::TripletList(int size) {
  ConstructTripletList_r_wrp(ih_this, &size);
}

////////////////////////////////////////////////////////////////////////////////
// TripletList::TripletList(const DistributedSparseMatrix &matrix) {
//   GetMatrixTripletList_wrp(&matrix.ih_this,&ih_this);
// }

////////////////////////////////////////////////////////////////////////////////
void TripletList::Resize(int size) { ResizeTripletList_r_wrp(ih_this, &size); }

////////////////////////////////////////////////////////////////////////////////
void TripletList::Append(const Triplet &value) {
  AppendToTripletList_r_wrp(ih_this, &(value.index_column), &(value.index_row),
                            &(value.point_value));
}

////////////////////////////////////////////////////////////////////////////////
void TripletList::SetTripletAt(int index, const Triplet &value) {
  int adjusted_index = index + 1;
  SetTripletAt_r_wrp(ih_this, &adjusted_index, &(value.index_column),
                     &(value.index_row), &(value.point_value));
}

////////////////////////////////////////////////////////////////////////////////
Triplet TripletList::GetTripletAt(int index) const {
  Triplet temp;
  int adjusted_index = index + 1;
  GetTripletAt_r_wrp(ih_this, &adjusted_index, &(temp.index_column),
                     &(temp.index_row), &(temp.point_value));
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
int TripletList::GetSize() const { return GetTripletListSize_r_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
TripletList::~TripletList() { DestructTripletList_r_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
void TripletList::SortTripletList(const TripletList &input, int matrix_columns,
                                  TripletList &sorted) {
  SortTripletList_r_wrp(input.ih_this, &matrix_columns, sorted.ih_this);
}
} // namespace NTPoly
