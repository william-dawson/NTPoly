#ifndef DISTRIBUTEDSPARSEMATRIX_h
#define DISTRIBUTEDSPARSEMATRIX_h

#include "Wrapper.h"
#include <string>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class DistributedMatrixMemoryPool;
class Permutation;
class SolverBase;
class TripletList_r;
class TripletList_c;
////////////////////////////////////////////////////////////////////////////////
//! A Module For Performing Distributed Sparse Matrix Operations.
class DistributedSparseMatrix {
public:
  //! Construct an empty matrix.
  //!\param matrix_dimension size fo the matrix.
  DistributedSparseMatrix(int matrix_dimension);
  //! Construct a matrix from file.
  //!\param file_name name of the file to build from.
  //!\param is_binary true if the file is a binary file.
  DistributedSparseMatrix(std::string file_name, bool is_binary = false);
  //! Copy constructor.
  //!\param matB matrix to opy from.
  DistributedSparseMatrix(const DistributedSparseMatrix &matB);

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

public:
  //! get the actual dimension of the matrix.
  int GetActualDimension() const;
  //! the logical dimension is scaled so each process has an even slice.
  int GetLogicalDimension() const;
  //! Extracts a triplet list of the data that is stored on this process.
  //! Data is returned with absolute coordinates.
  //! \param triplet_list the list to fill.
  void GetTripletList(TripletList_r &triplet_list);
  //! Extracts a triplet list of the data that is stored on this process.
  //! Data is returned with absolute coordinates.
  //! \param triplet_list the list to fill.
  void GetTripletList(TripletList_c &triplet_list);
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

public:
  //! Transpose a sparse matrix.
  //\param matA matrix to compute the transpose of.
  void Transpose(const NTPoly::DistributedSparseMatrix &matA);

public:
  //! this = dot(this,matB)
  //!\param matB the matrix to dot.
  double Dot(const DistributedSparseMatrix &matB);
  //! this = alpha*MatB + this (AXPY)
  //!\param matB the matrix to add.
  //!\param alpha scaling factor.
  //!\param threshold for flushing small values.
  void Increment(const DistributedSparseMatrix &matB, double alpha = 1.0,
                 double threshold = 0.0);
  //! Elementwise multiplication. C_ij = A_ij * B_ij. Also known as a Hadamard
  //  product.
  //!\param matA first mat.
  //!\param matB second mat.
  void PairwiseMultiply(const DistributedSparseMatrix &matA,
                        const DistributedSparseMatrix &matB);
  //! this := alpha*matA*matB+ beta*this (GEMM)
  //!\param matA first mat.
  //!\param matB second mat.
  //!\param memory_pool memory pool for intermediates.
  //!\param alpha scaling factor.
  //!\param beta scaling factor.
  //!\param threshold for flushing small values.
  void Gemm(const DistributedSparseMatrix &matA,
            const DistributedSparseMatrix &matB,
            DistributedMatrixMemoryPool &memory_pool, double alpha = 1.0,
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
  ~DistributedSparseMatrix();

private:
  int ih_this[SIZE_wrp];

private:
  //! Assignment operator, locked.
  DistributedSparseMatrix &operator=(const DistributedSparseMatrix &);
  friend class LoadBalancer;
  friend class SolverBase;
  template <class T> friend class TripletList;
  friend class DistributedMatrixMemoryPool;
};
} // namespace NTPoly
#endif
