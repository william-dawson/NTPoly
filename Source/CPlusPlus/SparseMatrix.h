#ifndef SparseMatrix_h
#define SparseMatrix_h

#include "Wrapper.h"
#include <string>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class MatrixMemoryPool;
class TripletList;
////////////////////////////////////////////////////////////////////////////////
//! A datatype for storing a CSR matrix.
class SparseMatrix {
public:
  //! Basic constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  SparseMatrix(int columns, int rows);
  //! Construct from a matrix market file.
  //!\param file_name matrix market file name.
  SparseMatrix(std::string file_name);
  //! Construct from a triplet list.
  //!\param list a list of triplet values to set in the matrix.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  SparseMatrix(const NTPoly::TripletList &list, int rows, int columns);
  //! Copy constructor.
  //!\param matB the matrix to copy from.
  SparseMatrix(const NTPoly::SparseMatrix &matB);

public:
  //! Get the number of rows in a matrix.
  int GetRows() const;
  //! Get the number of columns in a matrix.
  int GetColumns() const;

public:
  //! Scale the matrix by a constant.
  //!\param constant value to scale by.
  void Scale(double constant);
  //! This = alpha*MatrixB + This(AXPY).
  //!\param matB matrix to add.
  //!\param alpha scale for the matrix.
  //!\param threshold for flushing small values.
  void Increment(const NTPoly::SparseMatrix &matB, double alpha,
                 double threshold);
  //! Matrix dot product.
  //!\param matB matrix to dot with.
  //!\result the dot product of this and matB.
  double Dot(const NTPoly::SparseMatrix &matB) const;
  //! Pairwise multiply two sparse matrices.
  //!\param matA
  //!\param matB
  void PairwiseMultiply(const NTPoly::SparseMatrix &matA,
                        const NTPoly::SparseMatrix &matB);
  //! This := alpha*matA*op( matB ) + beta*this
  //!\param matA
  //!\param matB
  //!\param isATransposed true if A is already transposed.
  //!\param isBTransposed true if B is already transposed.
  //!\param alpha scaling value.
  //!\param beta scaling value.
  //!\param threshold for flushing small values.
  //!\param memory_pool a memory pool to use for storing intermediates.
  void Gemm(const NTPoly::SparseMatrix &matA, const NTPoly::SparseMatrix &matB,
            bool isATransposed, bool isBTransposed, double alpha, double beta,
            double threshold, NTPoly::MatrixMemoryPool &memory_pool);

public:
  //! Transpose a sparse matrix.
  //\param matA matrix to compute the transpose of.
  void Transpose(const NTPoly::SparseMatrix &matA);

public:
  //! Print the sparse matrix to the console.
  void Print();
  //! Write the sparse matrix to file.
  //!\param file_name file to print to.
  void WriteToMatrixMarket(std::string file_name);

public:
  //! Compute a triplet list from the entries in a matrix.
  //!\param triplet_list output.
  void MatrixToTripletList(NTPoly::TripletList &triplet_list);

public:
  //! Standard destructor.
  ~SparseMatrix();

private:
  //! Pointer to the underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  SparseMatrix &operator=(const SparseMatrix &);
};
} // namespace NTPoly
#endif
