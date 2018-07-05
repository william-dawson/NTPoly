#include "PMatrixMemoryPool.h"
#include "PSMatrix.h"
#include "Permutation.h"
#include "TripletList.h"
using std::string;
#include <iostream>
using std::cout;
using std::endl;
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "PSMatrix_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
Matrix_ps::Matrix_ps(int matrix_dimension) {
  ConstructEmptyMatrix_ps_wrp(ih_this, &matrix_dimension);
}

//////////////////////////////////////////////////////////////////////////////
Matrix_ps::Matrix_ps(std::string file_name, bool is_binary) {
  int string_length = file_name.length();
  if (is_binary) {
    ConstructMatrixFromBinary_ps_wrp(ih_this, &file_name.c_str()[0],
                                     &string_length);
  } else {
    ConstructMatrixFromMatrixMarket_ps_wrp(ih_this, &file_name.c_str()[0],
                                           &string_length);
  }
}

//////////////////////////////////////////////////////////////////////////////
Matrix_ps::Matrix_ps(const Matrix_ps &matB) {
  // Constructing empty here is important because the call to Empty_wrp
  // also allocates a handle.
  int matrix_dimension = matB.GetActualDimension();
  ConstructEmptyMatrix_ps_wrp(ih_this, &matrix_dimension);
  CopyMatrix_ps_wrp(matB.ih_this, ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::WriteToBinary(std::string file_name) const {
  int string_length = file_name.length();
  WriteMatrixToBinary_ps_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::WriteToMatrixMarket(string file_name) const {
  int string_length = file_name.length();
  WriteMatrixToMatrixMarket_ps_wrp(ih_this, &file_name.c_str()[0],
                                   &string_length);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::FillFromTripletList(const TripletList_r &triplet_list) {
  FillMatrixFromTripletList_psr_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::FillFromTripletList(const TripletList_c &triplet_list) {
  FillMatrixFromTripletList_psc_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::FillDistributedPermutation(const Permutation &lb,
                                           bool permuterows) {
  FillMatrixPermutation_ps_wrp(ih_this, lb.ih_this, &permuterows);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::FillIdentity() { FillMatrixIdentity_ps_wrp(ih_this); }

//////////////////////////////////////////////////////////////////////////////
int Matrix_ps::GetLogicalDimension() const {
  int temp;
  GetMatrixLogicalDimension_ps_wrp(ih_this, &temp);
  return temp;
}

//////////////////////////////////////////////////////////////////////////////
int Matrix_ps::GetActualDimension() const {
  int temp;
  GetMatrixActualDimension_ps_wrp(ih_this, &temp);
  return temp;
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::GetTripletList(TripletList_r &triplet_list) {
  GetMatrixTripletList_psr_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::GetTripletList(TripletList_c &triplet_list) {
  GetMatrixTripletList_psc_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::GetMatrixBlock(TripletList_r &triplet_list, int start_row,
                               int end_row, int start_column, int end_column) {
  GetMatrixBlock_psr_wrp(ih_this, triplet_list.ih_this, &start_row, &end_row,
                         &start_column, &end_column);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::GetMatrixBlock(TripletList_c &triplet_list, int start_row,
                               int end_row, int start_column, int end_column) {
  GetMatrixBlock_psc_wrp(ih_this, triplet_list.ih_this, &start_row, &end_row,
                         &start_column, &end_column);
}

////////////////////////////////////////////////////////////////////////////////
void Matrix_ps::Transpose(const Matrix_ps &matA) {
  TransposeMatrix_ps_wrp(matA.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void Matrix_ps::Conjugate() { ConjugateMatrix_ps_wrp(ih_this); }

//////////////////////////////////////////////////////////////////////////////
double Matrix_ps::Dot(const Matrix_ps &matB) {
  return DotMatrix_ps_wrp(ih_this, matB.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::Increment(const Matrix_ps &matB, double alpha,
                          double threshold) {
  IncrementMatrix_ps_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::PairwiseMultiply(const Matrix_ps &matA, const Matrix_ps &matB) {
  MatrixPairwiseMultiply_ps_wrp(matA.ih_this, matB.ih_this, ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::Gemm(const Matrix_ps &matA, const Matrix_ps &matB,
                     PMatrixMemoryPool &memory_pool, double alpha, double beta,
                     double threshold) {
  MatrixMultiply_ps_wrp(matA.ih_this, matB.ih_this, ih_this, &alpha, &beta,
                        &threshold, memory_pool.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void Matrix_ps::Scale(double constant) {
  ScaleMatrix_ps_wrp(ih_this, &constant);
}

//////////////////////////////////////////////////////////////////////////////
double Matrix_ps::Norm() const { return MatrixNorm_ps_wrp(ih_this); }

//////////////////////////////////////////////////////////////////////////////
double Matrix_ps::Trace() const {
  double temp_val;
  MatrixTrace_ps_wrp(ih_this, &temp_val);
  return temp_val;
}

//////////////////////////////////////////////////////////////////////////////
Matrix_ps::~Matrix_ps() { DestructMatrix_ps_wrp(ih_this); }
} // namespace NTPoly
