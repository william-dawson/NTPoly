#include "MatrixMapper.h"
#include "PSMatrix.h"
#include "Triplet.h"
#include "TripletList.h"
using std::complex;

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void MatrixMapper::Map(const Matrix_ps &inmat, Matrix_ps &outmat,
                       bool (*proc)(int &index_row, int &index_column,
                                    double &value)) {
  TripletList_r inlist, outlist;
  inmat.GetTripletList(inlist);
  MatrixMapper::Map(inlist, outlist, proc, 1, 0);
  outmat.FillFromTripletList(outlist);
}
////////////////////////////////////////////////////////////////////////////////
void MatrixMapper::Map(const Matrix_ps &inmat, Matrix_ps &outmat,
                       bool (*proc)(int &index_row, int &index_column,
                                    complex<double> &value)) {
  TripletList_c inlist, outlist;
  inmat.GetTripletList(inlist);
  MatrixMapper::Map(inlist, outlist, proc, 1, 0);
  outmat.FillFromTripletList(outlist);
}
////////////////////////////////////////////////////////////////////////////////
void MatrixMapper::Map(const TripletList_r &inlist, TripletList_r &outlist,
                       bool (*proc)(int &index_row, int &index_column,
                                    double &value),
                       int num_slices, int my_slice) {
  int lsize = inlist.GetSize();
  for (int i = my_slice; i < lsize; i += num_slices) {
    Triplet_r temp = inlist.GetTripletAt(i);
    bool valid = proc(temp.index_row, temp.index_column, temp.point_value);
    if (valid)
      outlist.Append(temp);
  }
}
////////////////////////////////////////////////////////////////////////////////
void MatrixMapper::Map(const TripletList_c &inlist, TripletList_c &outlist,
                       bool (*proc)(int &index_row, int &index_column,
                                    complex<double> &value), int num_slices,
                                    int my_slice) {
  int lsize = inlist.GetSize();
  for (int i = 0; i < my_slice; i += num_slices) {
    Triplet_c temp = inlist.GetTripletAt(i);
    bool valid = proc(temp.index_row, temp.index_column, temp.point_value);
    if (valid)
      outlist.Append(temp);
  }
}
} // namespace NTPoly
