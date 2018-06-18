#include "Triplet.h"
#include "TripletList.h"

#include <complex.h>

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "TripletList_c.h"
extern double _Complex real_to_complex(double x, double y);
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
TripletList_r::TripletList_r(int size) {
  ConstructTripletList_r_wrp(ih_this, &size);
}
TripletList_c::TripletList_c(int size) {
  ConstructTripletList_c_wrp(ih_this, &size);
}

////////////////////////////////////////////////////////////////////////////////
// TripletList::TripletList(const DistributedSparseMatrix &matrix) {
//   GetMatrixTripletList_wrp(&matrix.ih_this,&ih_this);
// }

////////////////////////////////////////////////////////////////////////////////
void TripletList_r::Resize(int size) {
  ResizeTripletList_r_wrp(ih_this, &size);
}
void TripletList_c::Resize(int size) {
  ResizeTripletList_c_wrp(ih_this, &size);
}

////////////////////////////////////////////////////////////////////////////////
void TripletList_r::Append(const Triplet_r &value) {
  AppendToTripletList_r_wrp(ih_this, &(value.index_column), &(value.index_row),
                            &(value.point_value));
}

void TripletList_c::Append(const Triplet_c &value) {
  double _Complex temp;
  temp = real_to_complex(value.point_value.real(), value.point_value.imag());
  AppendToTripletList_c_wrp(ih_this, &(value.index_column), &(value.index_row),
                            &temp);
}

////////////////////////////////////////////////////////////////////////////////

void TripletList_r::SetTripletAt(int index, const Triplet_r &value) {
  int adjusted_index = index + 1;
  SetTripletAt_r_wrp(ih_this, &adjusted_index, &(value.index_column),
                     &(value.index_row), &(value.point_value));
}

void TripletList_c::SetTripletAt(int index, const Triplet_c &value) {
  double _Complex temp;
  temp = real_to_complex(value.point_value.real(), value.point_value.imag());
  int adjusted_index = index + 1;
  SetTripletAt_c_wrp(ih_this, &adjusted_index, &(value.index_column),
                     &(value.index_row), &temp);
}

////////////////////////////////////////////////////////////////////////////////
Triplet_r TripletList_r::GetTripletAt(int index) const {
  Triplet_r temp;
  int adjusted_index = index + 1;
  GetTripletAt_r_wrp(ih_this, &adjusted_index, &(temp.index_column),
                     &(temp.index_row), &(temp.point_value));
  return temp;
}

Triplet_c TripletList_c::GetTripletAt(int index) const {
  Triplet_c temp;
  double _Complex tempc;
  int adjusted_index = index + 1;
  GetTripletAt_c_wrp(ih_this, &adjusted_index, &(temp.index_column),
                     &(temp.index_row), &tempc);
  temp.point_value = std::complex<double>(__real__ tempc, __imag__ tempc);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
int TripletList_r::GetSize() const { return GetTripletListSize_r_wrp(ih_this); }
int TripletList_c::GetSize() const { return GetTripletListSize_c_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
TripletList_r::~TripletList_r() { DestructTripletList_r_wrp(ih_this); }
TripletList_c::~TripletList_c() { DestructTripletList_c_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////

void TripletList_r::SortTripletList(const TripletList_r &input,
                                    int matrix_columns, TripletList_r &sorted) {
  SortTripletList_r_wrp(input.ih_this, &matrix_columns, sorted.ih_this);
}

void TripletList_c::SortTripletList(const TripletList_c &input,
                                    int matrix_columns, TripletList_c &sorted) {
  SortTripletList_c_wrp(input.ih_this, &matrix_columns, sorted.ih_this);
}
} // namespace NTPoly
