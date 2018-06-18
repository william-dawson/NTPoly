#ifndef SparseMatrix_h
#define SparseMatrix_h

#include "Wrapper.h"
#include <string>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class MatrixMemoryPool_r;
class TripletList_r;
class MatrixMemoryPool_c;
class TripletList_c;
////////////////////////////////////////////////////////////////////////////////
//! A datatype for storing a CSR matrix.
class SparseMatrix_r {
public:
  //! Basic constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  SparseMatrix_r(int columns, int rows);
  //! Construct from a matrix market file.
  //!\param file_name matrix market file name.
  SparseMatrix_r(std::string file_name);
  //! Construct from a triplet list.
  //!\param list a list of triplet values to set in the matrix.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  SparseMatrix_r(const NTPoly::TripletList_r &list, int rows, int columns);
  //! Copy constructor.
  //!\param matB the matrix to copy from.
  SparseMatrix_r(const NTPoly::SparseMatrix_r &matB);

public:
  //! Get the number of rows in a matrix.
  int GetRows() const;
  //! Get the number of columns in a matrix.
  int GetColumns() const;
  //! Extract a row from the matrix into the compressed vector representation.
  //!\param row_number the row to extract
  //!\param row_out the matrix representing that row
  void ExtractRow(int row_number, SparseMatrix_r &row_out) const;
  //! Extract a column from the matrix into the compressed vector representation
  //!\param column_number the column to extract
  //!\param column_out the matrix representing that column
  void ExtractColumn(int column_number, SparseMatrix_r &column_out) const;

public:
  //! Scale the matrix by a constant.
  //!\param constant value to scale by.
  void Scale(double constant);
  //! This = alpha*MatrixB + This(AXPY).
  //!\param matB matrix to add.
  //!\param alpha scale for the matrix.
  //!\param threshold for flushing small values.
  void Increment(const NTPoly::SparseMatrix_r &matB, double alpha,
                 double threshold);
  //! Matrix dot product.
  //!\param matB matrix to dot with.
  //!\result the dot product of this and matB.
  double Dot(const NTPoly::SparseMatrix_r &matB) const;
  //! Pairwise multiply two sparse matrices.
  //!\param matA
  //!\param matB
  void PairwiseMultiply(const NTPoly::SparseMatrix_r &matA,
                        const NTPoly::SparseMatrix_r &matB);
  //! This := alpha*matA*op( matB ) + beta*this
  //!\param matA
  //!\param matB
  //!\param isATransposed true if A is already transposed.
  //!\param isBTransposed true if B is already transposed.
  //!\param alpha scaling value.
  //!\param beta scaling value.
  //!\param threshold for flushing small values.
  //!\param memory_pool a memory pool to use for storing intermediates.
  void Gemm(const NTPoly::SparseMatrix_r &matA,
            const NTPoly::SparseMatrix_r &matB, bool isATransposed,
            bool isBTransposed, double alpha, double beta, double threshold,
            NTPoly::MatrixMemoryPool_r &memory_pool);

public:
  //! Compute the eigen vectors of a matrix.
  //!\param MatV the eigenvectors.
  //!\param threshold for pruning small values.
  void EigenDecomposition(NTPoly::SparseMatrix_r &MatV, double threshold) const;

public:
  //! Transpose a sparse matrix.
  //\param matA matrix to compute the transpose of.
  void Transpose(const NTPoly::SparseMatrix_r &matA);

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
  ~SparseMatrix_r();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  SparseMatrix &operator=(const SparseMatrix &);
};

class SparseMatrix_c {
public:
  //! Basic constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  SparseMatrix_c(int columns, int rows);
  //! Construct from a matrix market file.
  //!\param file_name matrix market file name.
  SparseMatrix_c(std::string file_name);
  //! Construct from a triplet list.
  //!\param list a list of triplet values to set in the matrix.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  SparseMatrix_c(const NTPoly::TripletList_c &list, int rows, int columns);
  //! Copy constructor.
  //!\param matB the matrix to copy from.
  SparseMatrix_c(const NTPoly::SparseMatrix_c &matB);

public:
  //! Get the number of rows in a matrix.
  int GetRows() const;
  //! Get the number of columns in a matrix.
  int GetColumns() const;
  //! Extract a row from the matrix into the compressed vector representation.
  //!\param row_number the row to extract
  //!\param row_out the matrix representing that row
  void ExtractRow(int row_number, SparseMatrix_c &row_out) const;
  //! Extract a column from the matrix into the compressed vector representation
  //!\param column_number the column to extract
  //!\param column_out the matrix representing that column
  void ExtractColumn(int column_number, SparseMatrix_c &column_out) const;

public:
  //! Scale the matrix by a constant.
  //!\param constant value to scale by.
  void Scale(double constant);
  //! This = alpha*MatrixB + This(AXPY).
  //!\param matB matrix to add.
  //!\param alpha scale for the matrix.
  //!\param threshold for flushing small values.
  void Increment(const NTPoly::SparseMatrix_c &matB, double alpha,
                 double threshold);
  //! Matrix dot product.
  //!\param matB matrix to dot with.
  //!\result the dot product of this and matB.
  double Dot(const NTPoly::SparseMatrix_c &matB) const;
  //! Pairwise multiply two sparse matrices.
  //!\param matA
  //!\param matB
  void PairwiseMultiply(const NTPoly::SparseMatrix_c &matA,
                        const NTPoly::SparseMatrix_c &matB);
  //! This := alpha*matA*op( matB ) + beta*this
  //!\param matA
  //!\param matB
  //!\param isATransposed true if A is already transposed.
  //!\param isBTransposed true if B is already transposed.
  //!\param alpha scaling value.
  //!\param beta scaling value.
  //!\param threshold for flushing small values.
  //!\param memory_pool a memory pool to use for storing intermediates.
  void Gemm(const NTPoly::SparseMatrix_c &matA,
            const NTPoly::SparseMatrix_c &matB, bool isATransposed,
            bool isBTransposed, double alpha, double beta, double threshold,
            NTPoly::MatrixMemoryPool_c &memory_pool);

public:
  //! Compute the eigen vectors of a matrix.
  //!\param MatV the eigenvectors.
  //!\param threshold for pruning small values.
  void EigenDecomposition(NTPoly::SparseMatrix_c &MatV, double threshold) const;

public:
  //! Transpose a sparse matrix.
  //\param matA matrix to compute the transpose of.
  void Transpose(const NTPoly::SparseMatrix_c &matA);
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
  ~SparseMatrix_c();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  SparseMatrix &operator=(const SparseMatrix &);
};
} // namespace NTPoly
#endif
