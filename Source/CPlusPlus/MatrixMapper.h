#ifndef MATRIXMAPS_h
#define MATRIXMAPS_h
#include <complex>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class Matrix_ps;
class TripletList_r;
class TripletList_c;
////////////////////////////////////////////////////////////////////////////////
class MatrixMapper {
public:
  //! Given a distributed matrix, apply this procedure to each element (real).
  //!\param inmat the matrix to apply the procedure to.
  //!\param outmat the matrix where each element has had proc called on it.
  //!\param proc the procedure to apply.
  static void Map(const Matrix_ps &inmat, Matrix_ps &outmat,
                  bool (*proc)(int &index_row, int &index_column,
                               double &value));
  //! Given a distributed matrix, apply this procedure to each element (cmplx).
  //!\param inmat the matrix to apply the procedure to.
  //!\param outmat the matrix where each element has had proc called on it.
  //!\param proc the procedure to apply.
  static void Map(const Matrix_ps &inmat, Matrix_ps &outmat,
                  bool (*proc)(int &index_row, int &index_column,
                               std::complex<double> &value));

private:
  //! Given a triplet list, apply this procedure to each element (real).
  //!\param inlist the list to apply the procedure to.
  //!\param outlist the list where each element has had proc called on it.
  //!\param proc the procedure to apply.
  static void Map(const TripletList_r &inlist, TripletList_r &outlist,
                  bool (*proc)(int &index_row, int &index_column,
                               double &value),
                  int num_slices = 1, int my_slice = 0);
  //! Given a triplet list, apply this procedure to each element (cmplx).
  //!\param inlist the list to apply the procedure to.
  //!\param outlist the list where each element has had proc called on it.
  //!\param proc the procedure to apply.
  static void Map(const TripletList_c &inmat, TripletList_c &outmat,
                  bool (*proc)(int &index_row, int &index_column,
                               std::complex<double> &value),
                  int num_slices = 1, int my_slice = 0);
  //! A helper that gets information about the process grid.
  //!\param mat the matrix to get the info of.
  //!\param num_slices how many slices is this matrix distributed on.
  //!\param my_slice which slice is this process on.
  static void GetSliceInfo(const Matrix_ps &mat, int &num_slices,
                           int &my_slice);
};
} // namespace NTPoly

#endif
