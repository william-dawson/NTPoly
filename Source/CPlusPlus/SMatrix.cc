#include "MatrixMemoryPool.h"
#include "SMatrix.h"
#include "TripletList.h"
using std::string;
using std::complex;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "SMatrix_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
Matrix_lsr::Matrix_lsr(int columns, int rows) {
  ConstructZeroMatrix_lsr_wrp(ih_this, &rows, &columns);
}
Matrix_lsc::Matrix_lsc(int columns, int rows) {
  ConstructZeroMatrix_lsc_wrp(ih_this, &rows, &columns);
}

////////////////////////////////////////////////////////////////////////////////
Matrix_lsr::Matrix_lsr(std::string file_name) {
  int string_length = file_name.length();
  ConstructMatrixFromFile_lsr_wrp(ih_this, &file_name.c_str()[0],
                                  &string_length);
}
Matrix_lsc::Matrix_lsc(std::string file_name) {
  int string_length = file_name.length();
  ConstructMatrixFromFile_lsc_wrp(ih_this, &file_name.c_str()[0],
                                  &string_length);
}

////////////////////////////////////////////////////////////////////////////////

Matrix_lsr::Matrix_lsr(const TripletList_r &list, int rows, int columns) {
  ConstructMatrixFromTripletList_lsr_wrp(ih_this, list.ih_this, &rows,
                                         &columns);
}

Matrix_lsc::Matrix_lsc(const TripletList_c &list, int rows, int columns) {
  ConstructMatrixFromTripletList_lsc_wrp(ih_this, list.ih_this, &rows,
                                         &columns);
}

////////////////////////////////////////////////////////////////////////////////

Matrix_lsr::Matrix_lsr(const Matrix_lsr &matB) {
  CopyMatrix_lsr_wrp(matB.ih_this, ih_this);
}

Matrix_lsc::Matrix_lsc(const Matrix_lsc &matB) {
  CopyMatrix_lsc_wrp(matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
Matrix_lsr::~Matrix_lsr() { DestructMatrix_lsr_wrp(ih_this); }
Matrix_lsc::~Matrix_lsc() { DestructMatrix_lsc_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int Matrix_lsr::GetRows() const {
  int temp;
  GetMatrixRows_lsr_wrp(ih_this, &temp);
  return temp;
}
int Matrix_lsc::GetRows() const {
  int temp;
  GetMatrixRows_lsc_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
int Matrix_lsr::GetColumns() const {
  int temp;
  GetMatrixColumns_lsr_wrp(ih_this, &temp);
  return temp;
}
int Matrix_lsc::GetColumns() const {
  int temp;
  GetMatrixColumns_lsc_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////

void Matrix_lsr::ExtractRow(int row_number, Matrix_lsr &row_out) const {
  int temp = row_number + 1;
  ExtractMatrixRow_lsr_wrp(ih_this, &temp, row_out.ih_this);
}

void Matrix_lsc::ExtractRow(int row_number, Matrix_lsc &row_out) const {
  int temp = row_number + 1;
  ExtractMatrixRow_lsc_wrp(ih_this, &temp, row_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////

void Matrix_lsr::ExtractColumn(int column_number, Matrix_lsr &column_out) const {
  int temp = column_number + 1;
  ExtractMatrixColumn_lsr_wrp(ih_this, &temp, column_out.ih_this);
}

void Matrix_lsc::ExtractColumn(int column_number, Matrix_lsc &column_out) const {
  int temp = column_number + 1;
  ExtractMatrixColumn_lsc_wrp(ih_this, &temp, column_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void Matrix_lsr::Scale(double constant) {
  ScaleMatrix_lsr_wrp(ih_this, &constant);
}
void Matrix_lsc::Scale(double constant) {
  ScaleMatrix_lsc_wrp(ih_this, &constant);
}

////////////////////////////////////////////////////////////////////////////////

void Matrix_lsr::Increment(const Matrix_lsr &matB, double alpha,
                          double threshold) {
  IncrementMatrix_lsr_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

void Matrix_lsc::Increment(const Matrix_lsc &matB, double alpha,
                          double threshold) {
  IncrementMatrix_lsc_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

////////////////////////////////////////////////////////////////////////////////

double Matrix_lsr::Dot(const Matrix_lsr &matB) const {
  double val;
  DotMatrix_lsr_wrp(ih_this, matB.ih_this, &val);
  return val;
}

complex<double> Matrix_lsc::Dot(const Matrix_lsc &matB) const {
  double real, imag;
  complex<double> val;
  DotMatrix_lsc_wrp(ih_this, matB.ih_this, &real, &imag);
  val = complex<double>(real,imag);
  return val;
}

////////////////////////////////////////////////////////////////////////////////

void Matrix_lsr::PairwiseMultiply(const NTPoly::Matrix_lsr &matA,
                                 const NTPoly::Matrix_lsr &matB) {
  PairwiseMultiplyMatrix_lsr_wrp(matA.ih_this, matB.ih_this, ih_this);
}

void Matrix_lsc::PairwiseMultiply(const NTPoly::Matrix_lsc &matA,
                                 const NTPoly::Matrix_lsc &matB) {
  PairwiseMultiplyMatrix_lsc_wrp(matA.ih_this, matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////

void Matrix_lsr::Gemm(const Matrix_lsr &matA, const Matrix_lsr &matB,
                     bool isATransposed, bool isBTransposed, double alpha,
                     double beta, double threshold,
                     MatrixMemoryPool_r &memory_pool) {
  MatrixMultiply_lsr_wrp(matA.ih_this, matB.ih_this, ih_this, &isATransposed,
                         &isBTransposed, &alpha, &beta, &threshold,
                         memory_pool.ih_this);
}

void Matrix_lsc::Gemm(const Matrix_lsc &matA, const Matrix_lsc &matB,
                     bool isATransposed, bool isBTransposed, double alpha,
                     double beta, double threshold,
                     MatrixMemoryPool_c &memory_pool) {
  MatrixMultiply_lsc_wrp(matA.ih_this, matB.ih_this, ih_this, &isATransposed,
                         &isBTransposed, &alpha, &beta, &threshold,
                         memory_pool.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void Matrix_lsr::Transpose(const Matrix_lsr &matA) {
  TransposeMatrix_lsr_wrp(matA.ih_this, ih_this);
}

void Matrix_lsc::Transpose(const Matrix_lsc &matA) {
  TransposeMatrix_lsc_wrp(matA.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void Matrix_lsc::Conjugate() { ConjugateMatrix_lsc_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
void Matrix_lsr::Print() const { PrintMatrix_lsr_wrp(ih_this); }
void Matrix_lsc::Print() const { PrintMatrix_lsc_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
void Matrix_lsr::WriteToMatrixMarket(string file_name) const {
  int string_length = file_name.length();
  PrintMatrixF_lsr_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

void Matrix_lsc::WriteToMatrixMarket(string file_name) const {
  int string_length = file_name.length();
  PrintMatrixF_lsc_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

////////////////////////////////////////////////////////////////////////////////

void Matrix_lsr::MatrixToTripletList(TripletList_r &triplet_list) const {
  MatrixToTripletList_lsr_wrp(ih_this, triplet_list.ih_this);
}

void Matrix_lsc::MatrixToTripletList(TripletList_c &triplet_list) const {
  MatrixToTripletList_lsc_wrp(ih_this, triplet_list.ih_this);
}

////////////////////////////////////////////////////////////////////////////////

void Matrix_lsr::EigenDecomposition(NTPoly::Matrix_lsr &MatV,
                                   double threshold) const {
  EigenDecomposition_lsr_wrp(ih_this, MatV.ih_this, &threshold);
}

void Matrix_lsc::EigenDecomposition(NTPoly::Matrix_lsc &MatV,
                                   double threshold) const {
  EigenDecomposition_lsc_wrp(ih_this, MatV.ih_this, &threshold);
}
} // namespace NTPoly
