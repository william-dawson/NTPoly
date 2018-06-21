!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE MatrixSModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixMarketModule, ONLY : ParseMMHeader
  USE MatrixModule, ONLY : Matrix_l
  USE VectorSModule, ONLY : AddSparseVectors, DotSparseVectors, &
      & PairwiseMultiplyVectors
  USE TripletListModule, ONLY : TripletList_r, TripletList_c
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, ABSTRACT, EXTENDS(Matrix_l), PUBLIC :: Matrix_ls
    INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
    INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
  END TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(Matrix_ls), PUBLIC :: Matrix_lsr
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values !< Values
  CONTAINS
    !! Construct/Destruct
    PROCEDURE :: InitEmpty => ConstructEmpty_lsr
    PROCEDURE :: InitFromFile => ConstructFromFile_lsr
    PROCEDURE :: InitFromTripletList => ConstructFromTripletList_lsr
    PROCEDURE :: Copy => Copy_lsr
    PROCEDURE :: Destruct => Destruct_lsr
    !! Basic Accessors
    PROCEDURE :: ExtractRow => ExtractRow_lsr
    PROCEDURE :: ExtractColumn => ExtractColumn_lsr
    !! Algebra
    PROCEDURE :: Scale => Scale_lsr
    PROCEDURE :: Increment => Increment_lsr
    PROCEDURE :: PairwiseMultiply => PairwiseMultiply_lsr
    PROCEDURE :: ColumnNorm => ColumnNorm_lsr
    PROCEDURE :: Norm => Norm_lsr
    !! ETC
    PROCEDURE :: Transpose => Transpose_lsr
    PROCEDURE :: Print => Print_lsr
    PROCEDURE :: PrintHeader => PrintHeader_lsr
  END TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(Matrix_ls), PUBLIC :: Matrix_lsc
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: values !< Values
  CONTAINS
    !! Construct/Destruct
    PROCEDURE :: InitEmpty => ConstructEmpty_lsc
    PROCEDURE :: InitFromFile => ConstructFromFile_lsc
    PROCEDURE :: InitFromTripletList => ConstructFromTripletList_lsc
    PROCEDURE :: Copy => Copy_lsc
    PROCEDURE :: Destruct => Destruct_lsc
    !! Basic Accessors
    PROCEDURE :: ExtractRow => ExtractRow_lsc
    PROCEDURE :: ExtractColumn => ExtractColumn_lsc
    !! Algebra
    PROCEDURE :: Scale => Scale_lsc
    PROCEDURE :: Increment => Increment_lsc
    PROCEDURE :: PairwiseMultiply => PairwiseMultiply_lsc
    PROCEDURE :: ColumnNorm => ColumnNorm_lsc
    PROCEDURE :: Norm => Norm_lsc
    !! ETC
    PROCEDURE :: Transpose => Transpose_lsc
    PROCEDURE :: Print => Print_lsc
    PROCEDURE :: PrintHeader => PrintHeader_lsc
  END TYPE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Internal only. Create a sparse matrix with a certain number of columns
  !! and rows. Will allocate storage for the outer values, nothing else.
  !! @param[out] this the matrix being created. It will have the outer
  !! index allocated, but nothing else.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE SUBROUTINE ConstructEmpty_lsr(this, rows, columns, zero_in)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: columns, rows
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    INCLUDE "sparse_includes/ConstructEmptyMatrix.f90"
  END SUBROUTINE ConstructEmpty_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Internal only. Create a sparse matrix with a certain number of columns
  !! and rows. Will allocate storage for the outer values, nothing else.
  !! @param[out] this the matrix being created. It will have the outer
  !! index allocated, but nothing else.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE SUBROUTINE ConstructEmpty_lsc(this, rows, columns, zero_in)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: columns, rows
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    INCLUDE "sparse_includes/ConstructEmptyMatrix.f90"
  END SUBROUTINE ConstructEmpty_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] this the matrix being constructed.
  !! @param[in] file_name name of the file.
  SUBROUTINE ConstructFromFile_lsr(this, file_name)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Triplet_r) :: temporary

     INCLUDE "sparse_includes/ConstructMatrixFromFile.f90"
  END SUBROUTINE ConstructFromFile_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] this the matrix being constructed.
  !! @param[in] file_name name of the file.
  SUBROUTINE ConstructFromFile_lsc(this, file_name)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Triplet_c) :: temporary

     INCLUDE "sparse_includes/ConstructMatrixFromFile.f90"
  END SUBROUTINE ConstructFromFile_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! The triplet list must be sorted to efficiently fill in the matrix. This
  !! constructor assumes \b you have already sorted the triplet list.
  !! @param[out] this the matrix being constructed
  !! @param[in] triplet_list a list of triplet values. They must be sorted.
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructFromTripletList_lsr(this, triplet_list, &
       & rows, columns)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(INOUT) :: this
    TYPE(TripletList_r), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN) :: rows, columns

    INCLUDE "sparse_includes/ConstructMatrixFromTripletList.f90"

  END SUBROUTINE ConstructFromTripletList_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! The triplet list must be sorted to efficiently fill in the matrix. This
  !! constructor assumes \b you have already sorted the triplet list.
  !! @param[out] this the matrix being constructed
  !! @param[in] triplet_list a list of triplet values. They must be sorted.
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructFromTripletList_lsc(this, triplet_list, &
       & rows, columns)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(INOUT) :: this
    TYPE(TripletList_c), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN) :: rows, columns

    INCLUDE "sparse_includes/ConstructMatrixFromTripletList.f90"

  END SUBROUTINE ConstructFromTripletList_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  !! @param[inout] this the matrix to free up
  PURE SUBROUTINE Destruct_lsr(this)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(INOUT) :: this

    INCLUDE "sparse_includes/DestructMatrix.f90"
  END SUBROUTINE Destruct_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  !! @param[inout] this the matrix to free up
  PURE SUBROUTINE Destruct_lsc(this)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(INOUT) :: this

    INCLUDE "sparse_includes/DestructMatrix.f90"
  END SUBROUTINE Destruct_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE Copy_lsr(matA, matB)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: matA
    CLASS(Matrix_l), INTENT(INOUT) :: matB

    SELECT TYPE(matB)
    CLASS IS(Matrix_lsr)
      INCLUDE "sparse_includes/CopyMatrix.f90"
    END SELECT
  END SUBROUTINE Copy_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE Copy_lsc(matA, matB)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: matA
    CLASS(Matrix_l), INTENT(INOUT) :: matB

    SELECT TYPE(matB)
    CLASS IS(Matrix_lsc)
      INCLUDE "sparse_includes/CopyMatrix.f90"
    END SELECT
  END SUBROUTINE Copy_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractRow_lsr(this, row_number, row_out)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    CLASS(Matrix_l), INTENT(INOUT) :: row_out
    !! Temporary Variables
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: value_buffer
    !! Temporary Variables
    INTEGER :: values_found
    INTEGER :: total_counter, elements_per_inner
    INTEGER :: outer_counter
    INTEGER :: inner_counter

    SELECT TYPE(row_out)
    CLASS IS(Matrix_lsr)
      INCLUDE "sparse_includes/ExtractMatrixRow.f90"
    END SELECT
  END SUBROUTINE ExtractRow_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractRow_lsc(this, row_number, row_out)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    CLASS(Matrix_l), INTENT(INOUT) :: row_out
    !! Temporary Variables
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: value_buffer
    !! Temporary Variables
    INTEGER :: values_found
    INTEGER :: total_counter, elements_per_inner
    INTEGER :: outer_counter
    INTEGER :: inner_counter

    SELECT TYPE(row_out)
    CLASS IS(Matrix_lsc)
      INCLUDE "sparse_includes/ExtractMatrixRow.f90"
    END SELECT
  END SUBROUTINE ExtractRow_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractColumn_lsr(this, column_number, column_out)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    CLASS(Matrix_l), INTENT(INOUT) :: column_out
    !! Local variables
    INTEGER :: number_of_values
    INTEGER :: start_index
    INTEGER :: counter

    SELECT TYPE(column_out)
    CLASS IS(Matrix_lsr)
      INCLUDE "sparse_includes/ExtractMatrixColumn.f90"
    END SELECT
  END SUBROUTINE ExtractColumn_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractColumn_lsc(this, column_number, column_out)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    CLASS(Matrix_l), INTENT(INOUT) :: column_out
    !! Local variables
    INTEGER :: number_of_values
    INTEGER :: start_index
    INTEGER :: counter

    SELECT TYPE(column_out)
    CLASS IS(Matrix_lsc)
      INCLUDE "sparse_includes/ExtractMatrixColumn.f90"
    END SELECT
  END SUBROUTINE ExtractColumn_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  !! @param[inout] matA Matrix A.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE Scale_lsr(matA,constant)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(INOUT) :: matA
    REAL(NTREAL), INTENT(IN) :: constant

    INCLUDE "sparse_includes/ScaleMatrix.f90"
  END SUBROUTINE Scale_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  !! @param[inout] matA Matrix A.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE Scale_lsc(matA,constant)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(INOUT) :: matA
    REAL(NTREAL), INTENT(IN) :: constant

    INCLUDE "sparse_includes/ScaleMatrix.f90"
  END SUBROUTINE Scale_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !! This will utilize the sparse vector addition routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier (optional, default=1.0)
  !! @param[in] threshold_in for flushing values to zero. (Optional, default=0).
  PURE SUBROUTINE Increment_lsr(matA, matB, alpha_in, threshold_in)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN)  :: matA
    CLASS(Matrix_l), INTENT(INOUT) :: matB
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_lsr) :: matC
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: size_of_a, size_of_b

    SELECT TYPE(matB)
    CLASS IS(Matrix_lsr)
      INCLUDE "sparse_includes/IncrementMatrix.f90"
    END SELECT
  END SUBROUTINE Increment_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !! This will utilize the sparse vector addition routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier (optional, default=1.0)
  !! @param[in] threshold_in for flushing values to zero. (Optional, default=0).
  PURE SUBROUTINE Increment_lsc(matA, matB, alpha_in, threshold_in)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN)  :: matA
    CLASS(Matrix_l), INTENT(INOUT) :: matB
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_lsc) :: matC
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: size_of_a, size_of_b

    SELECT TYPE(matB)
    CLASS IS(Matrix_lsc)
      INCLUDE "sparse_includes/IncrementMatrix.f90"
    END SELECT
  END SUBROUTINE Increment_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !! This will utilize the sparse vector pairwise routine.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] this = MatA mult MatB.
  PURE SUBROUTINE PairwiseMultiply_lsr(this, matA, matB)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(INOUT) :: this
    CLASS(Matrix_l), INTENT(IN)  :: matA
    CLASS(Matrix_l), INTENT(IN) :: matB
    !! Local Variables
    TYPE(Matrix_lsr) :: TempMat
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    INTEGER :: size_of_a, size_of_b

    SELECT TYPE(matA)
    CLASS IS(Matrix_lsr)
      SELECT TYPE(matB)
      CLASS IS(Matrix_lsr)
        INCLUDE "sparse_includes/PairwiseMultiplyMatrix.f90"
      END SELECT
    END SELECT
  END SUBROUTINE PairwiseMultiply_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !! This will utilize the sparse vector pairwise routine.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] this = MatA mult MatB.
  PURE SUBROUTINE PairwiseMultiply_lsc(this, matA, matB)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(INOUT) :: this
    CLASS(Matrix_l), INTENT(IN)  :: matA
    CLASS(Matrix_l), INTENT(IN) :: matB
    !! Local Variables
    TYPE(Matrix_lsc) :: TempMat
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    INTEGER :: size_of_a, size_of_b

    SELECT TYPE(matA)
    CLASS IS(Matrix_lsc)
      SELECT TYPE(matB)
      CLASS IS(Matrix_lsc)
        INCLUDE "sparse_includes/PairwiseMultiplyMatrix.f90"
      END SELECT
    END SELECT
  END SUBROUTINE PairwiseMultiply_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  !! @param[in] this the matrix to compute the norm of.
  !! @param[out] norm_per_column the norm value for each column in this matrix.
  PURE SUBROUTINE ColumnNorm_lsr(this, norm_per_column)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: this
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: norm_per_column
    !! Local Data
    REAL(NTREAL) :: temp_value

    INCLUDE "sparse_includes/MatrixColumnNorm.f90"
  END SUBROUTINE ColumnNorm_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  !! @param[in] this the matrix to compute the norm of.
  !! @param[out] norm_per_column the norm value for each column in this matrix.
  PURE SUBROUTINE ColumnNorm_lsc(this, norm_per_column)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: this
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: norm_per_column
    !! Local Data
    COMPLEX(NTCOMPLEX) :: temp_value

    INCLUDE "sparse_includes/MatrixColumnNorm.f90"
  END SUBROUTINE ColumnNorm_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a sparse matrix.
  !! @param[in] this the matrix to compute the norm of.
  !! @result norm the matrix.
  PURE FUNCTION Norm_lsr(this) RESULT(norm)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: this
    REAL(NTREAL) :: norm
    !! Local Variables
    REAL(NTREAL), DIMENSION(this%columns) :: column

    INCLUDE "sparse_includes/MatrixNorm.f90"

  END FUNCTION Norm_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a sparse matrix.
  !! @param[in] this the matrix to compute the norm of.
  !! @result norm the matrix.
  PURE FUNCTION Norm_lsc(this) RESULT(norm)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: this
    REAL(NTREAL) :: norm
    !! Local Variables
    REAL(NTREAL), DIMENSION(this%columns) :: column

    INCLUDE "sparse_includes/MatrixNorm.f90"

  END FUNCTION Norm_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! The current implementation has you go from matrix to triplet list,
  !! triplet list to transposed triplet list. The triplet list must then be
  !! sorted and then the return matrix is constructed.
  !! @param[in] this the matrix to be transposed.
  !! @param[out] matT the input matrix transposed.
  PURE SUBROUTINE Transpose_lsr(this, matT)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN)  :: this
    CLASS(Matrix_l), INTENT(INOUT) :: matT
    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    !! Temporary Variables
    INTEGER :: II, JJ
    INTEGER :: inner_index, insert_pt, this_offset
    INTEGER :: num_values, elements_per_inner

    SELECT TYPE(matT)
    CLASS IS (Matrix_lsr)
      INCLUDE "sparse_includes/TransposeMatrix.f90"
    END SELECT
  END SUBROUTINE Transpose_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! The current implementation has you go from matrix to triplet list,
  !! triplet list to transposed triplet list. The triplet list must then be
  !! sorted and then the return matrix is constructed.
  !! @param[in] this the matrix to be transposed.
  !! @param[out] matT the input matrix transposed.
  PURE SUBROUTINE Transpose_lsc(this, matT)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN)  :: this
    CLASS(Matrix_l), INTENT(INOUT) :: matT
    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    !! Temporary Variables
    INTEGER :: II, JJ
    INTEGER :: inner_index, insert_pt, this_offset
    INTEGER :: num_values, elements_per_inner

    SELECT TYPE(matT)
    CLASS IS (Matrix_lsc)
      INCLUDE "sparse_includes/TransposeMatrix.f90"
    END SELECT
  END SUBROUTINE Transpose_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE Print_lsr(this, file_name_in)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TripletList_r) :: triplet_list

    INCLUDE "sparse_includes/PrintMatrix.f90"

  END SUBROUTINE Print_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE Print_lsc(this, file_name_in)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TripletList_c) :: triplet_list

    INCLUDE "sparse_includes/PrintMatrix.f90"

  END SUBROUTINE Print_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the matrix market header for this matrix
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_handle an open file handle to write to.
  SUBROUTINE PrintHeader_lsr(this, file_handle)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    WRITE(file_handle,'(A)') "%%MatrixMarket matrix coordinate real general"

  END SUBROUTINE PrintHeader_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the matrix market header for this matrix
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_handle an open file handle to write to.
  SUBROUTINE PrintHeader_lsc(this, file_handle)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    WRITE(file_handle,'(A)') "%%MatrixMarket matrix coordinate complex general"

  END SUBROUTINE PrintHeader_lsc
END MODULE
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Create a big matrix from an array of matrices by putting them one next
!   !! to another.
!   !! @param[in] mat_array 2d array of matrices to compose.
!   !! @param[in] block_rows the number of rows of the array of blocks.
!   !! @param[in] block_columns the number of columns of the array of blocks.
!   !! @param[out] out_matrix the composed matrix.
!   PURE SUBROUTINE ComposeMatrix_lsr(mat_array, block_rows, block_columns, &
!        & out_matrix)
!     !! Parameters
!     CLASS(Matrix_lsr), DIMENSION(block_rows,block_columns), INTENT(IN) :: &
!          & mat_array
!     INTEGER, INTENT(IN) :: block_rows, block_columns
!     CLASS(Matrix_lsr), INTENT(INOUT) :: out_matrix
!     !! Local Data
!     CLASS(Matrix_lsr), DIMENSION(block_columns) :: merged_columns
!     CLASS(Matrix_lsr) :: Temp
!     CLASS(Matrix_lsr), DIMENSION(block_rows,block_columns) :: mat_t
!
!     INCLUDE "sparse_includes/ComposeMatrix.f90"
!   END SUBROUTINE ComposeMatrix_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
!   !! the first matrix are followed by the columns of the matrices in the list.
!   !! @param[in] mat_list list of matrices to compose.
!   !! @param[out] out_matrix = [Matrix 1 | Matrix 2, ...] .
!   PURE SUBROUTINE ComposeMatrixColumns_lsr(mat_list, out_matrix)
!     !! Parameters
!     CLASS(Matrix_lsr), DIMENSION(:), INTENT(IN) :: mat_list
!     CLASS(Matrix_lsr), INTENT(INOUT) :: out_matrix
!
!     INCLUDE "sparse_includes/ComposeMatrixColumns.f90"
!   END SUBROUTINE ComposeMatrixColumns_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Split a sparse matrix into an array of sparse matrices.
!   !! @param[in] this the matrix to split.
!   !! @param[in] block_rows number of rows to split the matrix into.
!   !! @param[in] block_columns number of columns to split the matrix into.
!   !! @param[out] split_array a COLUMNxROW array for the output to go into.
!   !! @param[in] block_size_row_in specifies the size of the  rows.
!   !! @param[in] block_size_column_in specifies the size of the columns.
!   SUBROUTINE SplitMatrix_lsr(this, block_rows, block_columns, &
!        & split_array, block_size_row_in, block_size_column_in)
!     !! Parameters
!     CLASS(Matrix_lsr), INTENT(IN) :: this
!     INTEGER, INTENT(IN) :: block_rows, block_columns
!     CLASS(Matrix_lsr), DIMENSION(:,:), INTENT(INOUT) :: split_array
!     INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
!     INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
!     !! Local Data
!     CLASS(Matrix_lsr), DIMENSION(block_columns) :: column_split
!     CLASS(Matrix_lsr), DIMENSION(block_rows) :: row_split
!     CLASS(Matrix_lsr) :: Temp
!
!     INCLUDE "sparse_includes/SplitMatrix.f90"
!   END SUBROUTINE SplitMatrix_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Split a matrix into into small blocks based on the specified offsets.
!   !! @param[in] this matrix to perform this operation on.
!   !! @param[in] num_blocks number of blocks to split into
!   !! @param[out] block_sizes the sizes used for splitting.
!   !! @param[out] split_list 1D array of blocks.
!   PURE SUBROUTINE SplitMatrixColumns_lsr(this, num_blocks, block_sizes, &
!        & split_list)
!     !! Parameters
!     CLASS(Matrix_lsr), INTENT(IN) :: this
!     INTEGER, INTENT(IN) :: num_blocks
!     INTEGER, DIMENSION(num_blocks), INTENT(IN) :: block_sizes
!     CLASS(Matrix_lsr), DIMENSION(num_blocks), INTENT(INOUT) :: split_list
!
!     INCLUDE "sparse_includes/SplitMatrixColumns.f90"
!   END SUBROUTINE SplitMatrixColumns_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Construct a triplet list from a matrix.
!   !! @param[in] this the matrix to construct the triplet list from.
!   !! @param[out] triplet_list the triplet list we created.
!   PURE SUBROUTINE MatrixToTripletList_lsr(this, triplet_list)
!     !! Parameters
!     CLASS(Matrix_lsr), INTENT(IN) :: this
!     TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
!     !! Local Variables
!     TYPE(Triplet_r) :: temporary
!
!     INCLUDE "sparse_includes/MatrixToTripletList.f90"
!   END SUBROUTINE MatrixToTripletList_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Every value in the matrix is changed into its complex conjugate.
!   !! @param[inout] this the matrix to compute the complex conjugate of.
!   PURE SUBROUTINE ConjugateMatrix_lsc(this)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(INOUT) :: this
!
!     this%values = CONJG(this%values)
!   END SUBROUTINE ConjugateMatrix_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END MODULE MatrixSModule
