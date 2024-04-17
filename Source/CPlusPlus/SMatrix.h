#ifndef SMatrix_h
#define SMatrix_h

#include "Wrapper.h"
#include <complex>
#include <string>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class MatrixMemoryPool_r;
class TripletList_r;
class MatrixMemoryPool_c;
class TripletList_c;
////////////////////////////////////////////////////////////////////////////////
//! A datatype for storing a CSR matrix (real values).
class Matrix_lsr {
public:
  //! Basic constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  Matrix_lsr(int columns, int rows);
  //! Construct from a matrix market file.
  //!\param file_name matrix market file name.
  Matrix_lsr(std::string file_name);
  //! Construct from a triplet list.
  //!\param list a list of triplet values to set in the matrix.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  Matrix_lsr(const NTPoly::TripletList_r &list, int rows, int columns);
  //! Copy constructor.
  //!\param matB the matrix to copy from.
  Matrix_lsr(const NTPoly::Matrix_lsr &matB);

public:
  //! Get the number of rows in a matrix.
  int GetRows() const;
  //! Get the number of columns in a matrix.
  int GetColumns() const;
  //! Extract a row from the matrix into the compressed vector representation.
  //!\param row_number the row to extract
  //!\param row_out the matrix representing that row
  void ExtractRow(int row_number, Matrix_lsr &row_out) const;
  //! Extract a column from the matrix into the compressed vector representation
  //!\param column_number the column to extract
  //!\param column_out the matrix representing that column
  void ExtractColumn(int column_number, Matrix_lsr &column_out) const;

public:
  //! Scale the matrix by a constant.
  //!\param constant value to scale by.
  void Scale(double constant);
  //! This = alpha*MatrixB + This(AXPY).
  //!\param matB matrix to add.
  //!\param alpha scale for the matrix.
  //!\param threshold for flushing small values.
  void Increment(const NTPoly::Matrix_lsr &matB, double alpha,
                 double threshold);
  //! Matrix dot product.
  //!\param matB matrix to dot with.
  //!\result the dot product of this and matB.
  double Dot(const NTPoly::Matrix_lsr &matB) const;
  //! Pairwise multiply two sparse matrices.
  //!\param matA
  //!\param matB
  void PairwiseMultiply(const NTPoly::Matrix_lsr &matA,
                        const NTPoly::Matrix_lsr &matB);
  //! This := alpha*matA*op( matB ) + beta*this
  //!\param matA
  //!\param matB
  //!\param isATransposed true if A is already transposed.
  //!\param isBTransposed true if B is already transposed.
  //!\param alpha scaling value.
  //!\param beta scaling value.
  //!\param threshold for flushing small values.
  //!\param memory_pool a memory pool to use for storing intermediates.
  void Gemm(const NTPoly::Matrix_lsr &matA, const NTPoly::Matrix_lsr &matB,
            bool isATransposed, bool isBTransposed, double alpha, double beta,
            double threshold, NTPoly::MatrixMemoryPool_r &memory_pool);
  //! Scale a matrix using a diagonal matrix (triplet list form).
  //!\param tlist the triplet list.
  //!\param threshold for flushing small values.
  void DiagonalScale(const NTPoly::TripletList_r &tlist);

public:
  //! Transpose a sparse matrix.
  //\param matA matrix to compute the transpose of.
  void Transpose(const NTPoly::Matrix_lsr &matA);

public:
  //! Print the sparse matrix to the console.
  void Print() const;
  //! Write the sparse matrix to file.
  //!\param file_name file to print to.
  void WriteToMatrixMarket(std::string file_name) const;

public:
  //! Compute a triplet list from the entries in a matrix.
  //!\param triplet_list output.
  void MatrixToTripletList(NTPoly::TripletList_r &triplet_list) const;

public:
  //! Standard destructor.
  ~Matrix_lsr();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  Matrix_lsr &operator=(const Matrix_lsr &);
};
//! A datatype for storing a CSR matrix (complex values).
class Matrix_lsc {
public:
  //! Basic constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  Matrix_lsc(int columns, int rows);
  //! Construct from a matrix market file.
  //!\param file_name matrix market file name.
  Matrix_lsc(std::string file_name);
  //! Construct from a triplet list.
  //!\param list a list of triplet values to set in the matrix.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  Matrix_lsc(const NTPoly::TripletList_c &list, int rows, int columns);
  //! Copy constructor.
  //!\param matB the matrix to copy from.
  Matrix_lsc(const NTPoly::Matrix_lsc &matB);

public:
  //! Get the number of rows in a matrix.
  int GetRows() const;
  //! Get the number of columns in a matrix.
  int GetColumns() const;
  //! Extract a row from the matrix into the compressed vector representation.
  //!\param row_number the row to extract
  //!\param row_out the matrix representing that row
  void ExtractRow(int row_number, Matrix_lsc &row_out) const;
  //! Extract a column from the matrix into the compressed vector representation
  //!\param column_number the column to extract
  //!\param column_out the matrix representing that column
  void ExtractColumn(int column_number, Matrix_lsc &column_out) const;

public:
  //! Scale the matrix by a constant.
  //!\param constant value to scale by.
  void Scale(double constant);
  //! This = alpha*MatrixB + This(AXPY).
  //!\param matB matrix to add.
  //!\param alpha scale for the matrix.
  //!\param threshold for flushing small values.
  void Increment(const NTPoly::Matrix_lsc &matB, double alpha,
                 double threshold);
  //! Matrix dot product.
  //!\param matB matrix to dot with.
  //!\result the dot product of this and matB.
  std::complex<double> Dot(const NTPoly::Matrix_lsc &matB) const;
  //! Pairwise multiply two sparse matrices.
  //!\param matA
  //!\param matB
  void PairwiseMultiply(const NTPoly::Matrix_lsc &matA,
                        const NTPoly::Matrix_lsc &matB);
  //! This := alpha*matA*op( matB ) + beta*this
  //!\param matA
  //!\param matB
  //!\param isATransposed true if A is already transposed.
  //!\param isBTransposed true if B is already transposed.
  //!\param alpha scaling value.
  //!\param beta scaling value.
  //!\param threshold for flushing small values.
  //!\param memory_pool a memory pool to use for storing intermediates.
  void Gemm(const NTPoly::Matrix_lsc &matA, const NTPoly::Matrix_lsc &matB,
            bool isATransposed, bool isBTransposed, double alpha, double beta,
            double threshold, NTPoly::MatrixMemoryPool_c &memory_pool);
  //! Scale a matrix using a diagonal matrix (triplet list form).
  //!\param tlist the triplet list.
  //!\param threshold for flushing small values.
  void DiagonalScale(const NTPoly::TripletList_c &tlist, double threshold);

public:
  //! Transpose a sparse matrix.
  //\param matA matrix to compute the transpose of.
  void Transpose(const NTPoly::Matrix_lsc &matA);
  //! Compute the complex conjugate of a matrix
  void Conjugate();

public:
  //! Print the sparse matrix to the console.
  void Print() const;
  //! Write the sparse matrix to file.
  //!\param file_name file to print to.
  void WriteToMatrixMarket(std::string file_name) const;

public:
  //! Compute a triplet list from the entries in a matrix.
  //!\param triplet_list output.
  void MatrixToTripletList(NTPoly::TripletList_c &triplet_list) const;

public:
  //! Standard destructor.
  ~Matrix_lsc();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  Matrix_lsc &operator=(const Matrix_lsc &);
};
} // namespace NTPoly
#endif
