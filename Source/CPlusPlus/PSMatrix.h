#ifndef DISTRIBUTEDSPARSEMATRIX_h
#define DISTRIBUTEDSPARSEMATRIX_h

#include "Wrapper.h"
#include <complex>
#include <string>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class PMatrixMemoryPool;
class Permutation;
class SolverBase;
class TripletList_r;
class TripletList_c;
class ProcessGrid;
class MatrixMapper;
class MatrixConversion;
////////////////////////////////////////////////////////////////////////////////
//! A Module For Performing Distributed Sparse Matrix Operations.
class Matrix_ps {
public:
  //! Construct an empty matrix.
  //!\param matrix_dimension size fo the matrix.
  Matrix_ps(int matrix_dimension);
  //! Construct an empty matrix.
  //!\param matrix_dimension size fo the matrix.
  //!\param grid the process grid this matrix is distributed on.
  Matrix_ps(int matrix_dimension, const ProcessGrid &grid);
  //! Construct a matrix from file.
  //!\param file_name name of the file to build from.
  //!\param is_binary true if the file is a binary file.
  Matrix_ps(std::string file_name, bool is_binary = false);
  //! Construct a matrix from file.
  //!\param file_name name of the file to build from.
  //!\param is_binary true if the file is a binary file.
  //!\param grid the process grid this matrix is distributed on.
  Matrix_ps(std::string file_name, const ProcessGrid &grid,
            bool is_binary = false);
  //! Copy constructor.
  //!\param matB matrix to copy from.
  Matrix_ps(const Matrix_ps &matB);

public:
  //! Write the matrix to a custom binary format.
  //!\param file_name file to write to.
  void WriteToBinary(std::string file_name) const;
  //! Write the matrix to the matrix market format.
  //!\param file_name file to write to.
  void WriteToMatrixMarket(std::string file_name) const;

public:
  //! Fill in the matrix based on the contents of triplet lists.
  //!\param triplet_list list of values. Need to be absolute coordinates.
  void FillFromTripletList(const TripletList_r &triplet_list);
  //! Fill in the matrix based on the contents of triplet lists.
  //!\param triplet_list list of values. Need to be absolute coordinates.
  void FillFromTripletList(const TripletList_c &triplet_list);
  //! Fill the matrix based on a permutation.
  //!\param lb the permutation.
  //!\param permuterows true if this is a row permutation matrix.
  void FillDistributedPermutation(const Permutation &lb,
                                  bool permuterows = true);
  //! Fills this matrix as the identity matrix.
  void FillIdentity();
  //! This routine will fill a dense matrix so that every element has a given
  //! a value of 1. This is useful as a starting point for further filtering
  //! or mapping operations.
  void FillDense();

public:
  //! get the actual dimension of the matrix.
  int GetActualDimension() const;
  //! the logical dimension is scaled so each process has an even slice.
  int GetLogicalDimension() const;
  //! Get the total number of non-zero entries in the matrix
  long int GetSize() const;
  //! Extracts a triplet list of the data that is stored on this process.
  //! Data is returned with absolute coordinates.
  //! \param triplet_list the list to fill.
  void GetTripletList(TripletList_r &triplet_list) const;
  //! Extracts a triplet list of the data that is stored on this process.
  //! Data is returned with absolute coordinates.
  //! \param triplet_list the list to fill.
  void GetTripletList(TripletList_c &triplet_list) const;
  //! Extract an arbitrary block of a matrix into a triplet list. Block is
  //! defined by the row/column start/end values.
  //! This is slower than GetTripletList, because communication is required.
  //! Data is returned with absolute coordinates.
  //! \param triplet_list the list to fill.
  //! \param start_row the starting row for data to store on this process.
  //! \param end_row the ending row for data to store on this process.
  //! \param start_column the starting col for data to store on this process
  //! \param end_column the ending col for data to store on this process
  void GetMatrixBlock(TripletList_r &triplet_list, int start_row, int end_row,
                      int start_column, int end_column);
  //! Extract an arbitrary block of a matrix into a triplet list. Block is
  //! defined by the row/column start/end values.
  //! This is slower than GetTripletList, because communication is required.
  //! Data is returned with absolute coordinates.
  //! \param triplet_list the list to fill.
  //! \param start_row the starting row for data to store on this process.
  //! \param end_row the ending row for data to store on this process.
  //! \param start_column the starting col for data to store on this process
  //! \param end_column the ending col for data to store on this process
  void GetMatrixBlock(TripletList_c &triplet_list, int start_row, int end_row,
                      int start_column, int end_column);
  //! Copy an arbitrary slice from a matrix into a new smaller matrix.
  //! NTPoly only works with square matrices, so if the number of rows and
  //! columns is different the matrix is resized to the maximum size.
  //! \param submatrix the slice to fill.
  //! \param start_row the starting row to include in this matrix.
  //! \param end_row the ending row to include in this matrix.
  //! \param start_column the starting column to include in this matrix.
  //! \param end_column the last column to include in this matrix.
  void GetMatrixSlice(Matrix_ps &submatrix, int start_row, int end_row,
                      int start_column, int end_column);

public:
  //! Transpose a sparse matrix.
  //\param matA matrix to compute the transpose of.
  void Transpose(const NTPoly::Matrix_ps &matA);
  //! Compute the complex conjugate of a matrix
  void Conjugate();

public:
  //! Change the size of a matrix.
  //! If the new size is smaller, then values outside that range are deleted.
  //! IF the new size is bigger, zero padding is applied.
  //! Warning: this requires a full data redistribution.
  //\param new_size the new size of the matrix.
  void Resize(int new_size);

public:
  //! this = dot(this,matB)
  //!\param matB the matrix to dot.
  double Dot(const Matrix_ps &matB);
  //! this = dot(this,matB)
  //!\param matB the matrix to dot.
  std::complex<double> Dot_c(const Matrix_ps &matB);
  //! this = alpha*MatB + this (AXPY)
  //!\param matB the matrix to add.
  //!\param alpha scaling factor.
  //!\param threshold for flushing small values.
  void Increment(const Matrix_ps &matB, double alpha = 1.0,
                 double threshold = 0.0);
  //! Elementwise multiplication. C_ij = A_ij * B_ij. Also known as a Hadamard
  //  product.
  //!\param matA first mat.
  //!\param matB second mat.
  void PairwiseMultiply(const Matrix_ps &matA, const Matrix_ps &matB);
  //! this := alpha*matA*matB+ beta*this (GEMM)
  //!\param matA first mat.
  //!\param matB second mat.
  //!\param memory_pool memory pool for intermediates.
  //!\param alpha scaling factor.
  //!\param beta scaling factor.
  //!\param threshold for flushing small values.
  void Gemm(const Matrix_ps &matA, const Matrix_ps &matB,
            PMatrixMemoryPool &memory_pool, double alpha = 1.0,
            double beta = 0.0, double threshold = 0.0);
  //! scale the matrix by a constatn.
  //! constant the value to scale by.
  void Scale(double constant);
  //! compute the norm of a matrix.
  double Norm() const;
  //! compute the trace of a matrix.
  double Trace() const;

public:
  //! Destructor.
  ~Matrix_ps();

private:
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  Matrix_ps &operator=(const Matrix_ps &);
  friend class LoadBalancer;
  friend class SolverBase;
  template <class T> friend class TripletList;
  friend class PMatrixMemoryPool;
  friend class MatrixMapper;
  friend class MatrixConversion;
};
} // namespace NTPoly
#endif
