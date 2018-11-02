#include "MatrixMapper.h"
#include "PSMatrix.h"
#include "Triplet.h"
#include "TripletList.h"
#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "PSMatrix_c.h"
#include "ProcessGrid_c.h"
}

using std::complex;
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void MatrixMapper::GetSliceInfo(const Matrix_ps &mat, int &num_slices,
                                int &my_slice) {
  int ih_grid[SIZE_wrp];
  GetMatrixProcessGrid_ps_wrp(mat.ih_this, ih_grid);
  num_slices = GetNumSlices_wrp(ih_grid);
  my_slice = GetMySlice_wrp(ih_grid);
}
////////////////////////////////////////////////////////////////////////////////
void MatrixMapper::Map(const Matrix_ps &inmat, Matrix_ps &outmat,
                       bool (*proc)(int &index_row, int &index_column,
                                    double &value)) {
  int num_slices, my_slice;
  GetSliceInfo(inmat, num_slices, my_slice);
  TripletList_r inlist, outlist;
  inmat.GetTripletList(inlist);
  MatrixMapper::Map(inlist, outlist, proc, num_slices, my_slice);
  outmat.FillFromTripletList(outlist);
}
////////////////////////////////////////////////////////////////////////////////
void MatrixMapper::Map(const Matrix_ps &inmat, Matrix_ps &outmat,
                       bool (*proc)(int &index_row, int &index_column,
                                    complex<double> &value)) {
  int num_slices, my_slice;
  GetSliceInfo(inmat, num_slices, my_slice);
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
                                    complex<double> &value),
                       int num_slices, int my_slice) {
  int lsize = inlist.GetSize();
  for (int i = 0; i < my_slice; i += num_slices) {
    Triplet_c temp = inlist.GetTripletAt(i);
    bool valid = proc(temp.index_row, temp.index_column, temp.point_value);
    if (valid)
      outlist.Append(temp);
  }
}
} // namespace NTPoly
