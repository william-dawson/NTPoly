!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE MatrixSModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE MatrixMarketModule, ONLY : ParseMMHeader
  USE TripletListModule, ONLY: TripletList_r, TripletList_c, SortTripletList, &
       & DestructTripletList, SetTripletAt, ConstructTripletList, &
       & AppendToTripletList, SymmetrizeTripletList, ConvertTripletListType
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a local, real CSR matrix.
  TYPE, PUBLIC :: Matrix_lsr
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a local, complex CSR matrix.
  TYPE, PUBLIC :: Matrix_lsc
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Construct/Destruct
  PUBLIC :: ConstructEmptyMatrix
  PUBLIC :: ConstructMatrixFromTripletList
  PUBLIC :: DestructMatrix
  PUBLIC :: CopyMatrix
  !! Basic Accessors
  PUBLIC :: GetMatrixRows
  PUBLIC :: GetMatrixColumns
  PUBLIC :: ExtractMatrixRow
  PUBLIC :: ExtractMatrixColumn
  !! Routines for splitting and composing
  PUBLIC :: SplitMatrix
  PUBLIC :: SplitMatrixColumns
  PUBLIC :: ComposeMatrix
  PUBLIC :: ComposeMatrixColumns
  !! ETC
  PUBLIC :: ConvertMatrixType
  PUBLIC :: TransposeMatrix
  PUBLIC :: ConjugateMatrix
  PUBLIC :: PrintMatrix
  PUBLIC :: MatrixToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE Matrix_lsr
     MODULE PROCEDURE ConstructEmptyMatrix_lsr
     MODULE PROCEDURE ConstructMatrixFromFile_lsr
     MODULE PROCEDURE ConstructMatrixFromTripletList_lsr
  END INTERFACE
  INTERFACE Matrix_lsc
     MODULE PROCEDURE ConstructEmptyMatrix_lsc
     MODULE PROCEDURE ConstructMatrixFromFile_lsc
     MODULE PROCEDURE ConstructMatrixFromTripletList_lsc
  END INTERFACE
  INTERFACE ConstructEmptyMatrix
     MODULE PROCEDURE ConstructEmptyMatrixSub_lsr
     MODULE PROCEDURE ConstructEmptyMatrixSub_lsc
  END INTERFACE
  INTERFACE ConstructMatrixFromFile
     MODULE PROCEDURE ConstructMatrixFromFileSub_lsr
     MODULE PROCEDURE ConstructMatrixFromFileSub_lsc
  END INTERFACE
  INTERFACE ConstructMatrixFromTripletList
     MODULE PROCEDURE ConstructMatrixFromTripletListSub_lsr
     MODULE PROCEDURE ConstructMatrixFromTripletListSub_lsc
  END INTERFACE
  INTERFACE DestructMatrix
     MODULE PROCEDURE DestructMatrix_lsr
     MODULE PROCEDURE DestructMatrix_lsc
  END INTERFACE
  INTERFACE CopyMatrix
     MODULE PROCEDURE CopyMatrix_lsr
     MODULE PROCEDURE CopyMatrix_lsc
  END INTERFACE
  INTERFACE GetMatrixRows
     MODULE PROCEDURE GetMatrixRows_lsr
     MODULE PROCEDURE GetMatrixRows_lsc
  END INTERFACE
  INTERFACE GetMatrixColumns
     MODULE PROCEDURE GetMatrixColumns_lsr
     MODULE PROCEDURE GetMatrixColumns_lsc
  END INTERFACE
  INTERFACE ExtractMatrixRow
     MODULE PROCEDURE ExtractMatrixRow_lsr
     MODULE PROCEDURE ExtractMatrixRow_lsc
  END INTERFACE
  INTERFACE ExtractMatrixColumn
     MODULE PROCEDURE ExtractMatrixColumn_lsr
     MODULE PROCEDURE ExtractMatrixColumn_lsc
  END INTERFACE
  INTERFACE SplitMatrix
     MODULE PROCEDURE SplitMatrix_lsr
     MODULE PROCEDURE SplitMatrix_lsc
  END INTERFACE
  INTERFACE SplitMatrixColumns
     MODULE PROCEDURE SplitMatrixColumns_lsr
     MODULE PROCEDURE SplitMatrixColumns_lsc
  END INTERFACE
  INTERFACE ComposeMatrix
     MODULE PROCEDURE ComposeMatrix_lsr
     MODULE PROCEDURE ComposeMatrix_lsc
  END INTERFACE
  INTERFACE ComposeMatrixColumns
     MODULE PROCEDURE ComposeMatrixColumns_lsr
     MODULE PROCEDURE ComposeMatrixColumns_lsc
  END INTERFACE
  INTERFACE TransposeMatrix
     MODULE PROCEDURE TransposeMatrix_lsr
     MODULE PROCEDURE TransposeMatrix_lsc
  END INTERFACE
  INTERFACE ConjugateMatrix
     MODULE PROCEDURE ConjugateMatrix_lsc
  END INTERFACE
  INTERFACE PrintMatrix
     MODULE PROCEDURE PrintMatrix_lsr
     MODULE PROCEDURE PrintMatrix_lsc
  END INTERFACE
  INTERFACE MatrixToTripletList
     MODULE PROCEDURE MatrixToTripletList_lsr
     MODULE PROCEDURE MatrixToTripletList_lsc
  END INTERFACE
  INTERFACE ConvertMatrixType
     MODULE PROCEDURE ConvertMatrixType_lsrtolsc
     MODULE PROCEDURE ConvertMatrixType_lsctolsr
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE ConstructEmptyMatrixSub_lsr(this, rows, columns, zero_in)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: columns, rows
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    IF (PRESENT(zero_in)) THEN
       this = ConstructEmptyMatrix_lsr(rows, columns, zero_in)
    ELSE
       this = ConstructEmptyMatrix_lsr(rows, columns)
    ENDIF

  END SUBROUTINE ConstructEmptyMatrixSub_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Internal only. Create a sparse matrix with a certain number of columns
  !! and rows. Will allocate storage for the outer values, nothing else.
  !! @param[out] this the matrix being created. It will have the outer
  !! index allocated, but nothing else.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE FUNCTION ConstructEmptyMatrix_lsr(rows, columns, zero_in) RESULT(this)
    !! Parameters
    TYPE(Matrix_lsr) :: this
    INTEGER, INTENT(IN) :: columns, rows
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    INCLUDE "sparse_includes/ConstructEmptyMatrix.f90"
  END FUNCTION ConstructEmptyMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ConstructMatrixFromFileSub_lsr(this, file_name)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name

    this = ConstructMatrixFromFile_lsr(file_name)
  END SUBROUTINE ConstructMatrixFromFileSub_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] this the matrix being constructed.
  !! @param[in] file_name name of the file.
  FUNCTION ConstructMatrixFromFile_lsr(file_name) RESULT(this)
    !! Parameters
    TYPE(Matrix_lsr) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Triplet_r) :: temporary

#include "sparse_includes/ConstructMatrixFromFile.f90"
  END FUNCTION ConstructMatrixFromFile_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE ConstructMatrixFromTripletListSub_lsr(this, triplet_list, &
       & rows, columns)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(INOUT) :: this
    TYPE(TripletList_r), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN) :: rows, columns

    this = ConstructMatrixFromTripletList_lsr(triplet_list, rows, columns)

  END SUBROUTINE ConstructMatrixFromTripletListSub_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! The triplet list must be sorted to efficiently fill in the matrix. This
  !! constructor assumes \b you have already sorted the triplet list.
  !! @param[out] this the matrix being constructed
  !! @param[in] triplet_list a list of triplet values. They must be sorted.
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE FUNCTION ConstructMatrixFromTripletList_lsr(triplet_list,rows,columns) &
       & RESULT(this)
    !! Parameters
    TYPE(Matrix_lsr) :: this
    TYPE(TripletList_r), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN) :: rows, columns

#define ISCOMPLEX
#include "sparse_includes/ConstructMatrixFromTripletList.f90"
#undef ISCOMPLEX

  END FUNCTION ConstructMatrixFromTripletList_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  !! @param[inout] this the matrix to free up
  PURE SUBROUTINE DestructMatrix_lsr(this)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(INOUT) :: this

    INCLUDE "sparse_includes/DestructMatrix.f90"
  END SUBROUTINE DestructMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE CopyMatrix_lsr(matA, matB)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: matA
    TYPE(Matrix_lsr), INTENT(INOUT) :: matB

    INCLUDE "sparse_includes/CopyMatrix.f90"
  END SUBROUTINE CopyMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  !! @param[in] this the matrix.
  !! @result number of rows.
  PURE FUNCTION GetMatrixRows_lsr(this) RESULT(rows)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this

    INCLUDE "sparse_includes/GetMatrixRows.f90"
  END FUNCTION GetMatrixRows_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  !! @param[in] this the matrix.
  !! @result number of columns.
  PURE FUNCTION GetMatrixColumns_lsr(this) RESULT(columns)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this

    INCLUDE "sparse_includes/GetMatrixColumns.f90"
  END FUNCTION GetMatrixColumns_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixRow_lsr(this, row_number, row_out)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    TYPE(Matrix_lsr), INTENT(INOUT) :: row_out
    !! Temporary Variables
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: value_buffer

    INCLUDE "sparse_includes/ExtractMatrixRow.f90"
  END SUBROUTINE ExtractMatrixRow_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixColumn_lsr(this, column_number, column_out)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    TYPE(Matrix_lsr), INTENT(INOUT) :: column_out

    INCLUDE "sparse_includes/ExtractMatrixColumn.f90"
  END SUBROUTINE ExtractMatrixColumn_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! The current implementation has you go from matrix to triplet list,
  !! triplet list to transposed triplet list. The triplet list must then be
  !! sorted and then the return matrix is constructed.
  !! @param[in] this the matrix to be transposed.
  !! @param[out] matT the input matrix transposed.
  PURE SUBROUTINE TransposeMatrix_lsr(this, matT)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)  :: this
    TYPE(Matrix_lsr), INTENT(INOUT) :: matT

    INCLUDE "sparse_includes/TransposeMatrix.f90"
  END SUBROUTINE TransposeMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !! to another.
  !! @param[in] mat_array 2d array of matrices to compose.
  !! @param[in] block_rows the number of rows of the array of blocks.
  !! @param[in] block_columns the number of columns of the array of blocks.
  !! @param[out] out_matrix the composed matrix.
  PURE SUBROUTINE ComposeMatrix_lsr(mat_array, block_rows, block_columns, &
       & out_matrix)
    !! Parameters
    TYPE(Matrix_lsr), DIMENSION(block_rows,block_columns), INTENT(IN) :: &
         & mat_array
    INTEGER, INTENT(IN) :: block_rows, block_columns
    TYPE(Matrix_lsr), INTENT(INOUT) :: out_matrix
    !! Local Data
    TYPE(Matrix_lsr), DIMENSION(block_columns) :: merged_columns
    TYPE(Matrix_lsr) :: Temp
    TYPE(Matrix_lsr), DIMENSION(block_rows,block_columns) :: mat_t

    INCLUDE "sparse_includes/ComposeMatrix.f90"
  END SUBROUTINE ComposeMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
  !! the first matrix are followed by the columns of the matrices in the list.
  !! @param[in] mat_list list of matrices to compose.
  !! @param[out] out_matrix = [Matrix 1 | Matrix 2, ...] .
  PURE SUBROUTINE ComposeMatrixColumns_lsr(mat_list, out_matrix)
    !! Parameters
    TYPE(Matrix_lsr), DIMENSION(:), INTENT(IN) :: mat_list
    TYPE(Matrix_lsr), INTENT(INOUT) :: out_matrix

    INCLUDE "sparse_includes/ComposeMatrixColumns.f90"
  END SUBROUTINE ComposeMatrixColumns_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  !! @param[in] this the matrix to split.
  !! @param[in] block_rows number of rows to split the matrix into.
  !! @param[in] block_columns number of columns to split the matrix into.
  !! @param[out] split_array a COLUMNxROW array for the output to go into.
  !! @param[in] block_size_row_in specifies the size of the  rows.
  !! @param[in] block_size_column_in specifies the size of the columns.
  SUBROUTINE SplitMatrix_lsr(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: block_rows, block_columns
    TYPE(Matrix_lsr), DIMENSION(:,:), INTENT(INOUT) :: split_array
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
    !! Local Data
    TYPE(Matrix_lsr), DIMENSION(block_columns) :: column_split
    TYPE(Matrix_lsr), DIMENSION(block_rows) :: row_split
    TYPE(Matrix_lsr) :: Temp

    INCLUDE "sparse_includes/SplitMatrix.f90"
  END SUBROUTINE SplitMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a matrix into into small blocks based on the specified offsets.
  !! @param[in] this matrix to perform this operation on.
  !! @param[in] num_blocks number of blocks to split into
  !! @param[out] block_sizes the sizes used for splitting.
  !! @param[out] split_list 1D array of blocks.
  PURE SUBROUTINE SplitMatrixColumns_lsr(this, num_blocks, block_sizes, &
       & split_list)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: num_blocks
    INTEGER, DIMENSION(num_blocks), INTENT(IN) :: block_sizes
    TYPE(Matrix_lsr), DIMENSION(num_blocks), INTENT(INOUT) :: split_list

    INCLUDE "sparse_includes/SplitMatrixColumns.f90"
  END SUBROUTINE SplitMatrixColumns_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list from a matrix.
  !! @param[in] this the matrix to construct the triplet list from.
  !! @param[out] triplet_list the triplet list we created.
  PURE SUBROUTINE MatrixToTripletList_lsr(this, triplet_list)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    !! Local Variables
    TYPE(Triplet_r) :: temporary

    INCLUDE "sparse_includes/MatrixToTripletList.f90"
  END SUBROUTINE MatrixToTripletList_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintMatrix_lsr(this, file_name_in)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TripletList_r) :: triplet_list

#include "sparse_includes/PrintMatrix.f90"

  END SUBROUTINE PrintMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE ConstructEmptyMatrixSub_lsc(this, rows, columns, zero_in)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: columns, rows
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    IF (PRESENT(zero_in)) THEN
       this = ConstructEmptyMatrix_lsc(rows, columns, zero_in)
    ELSE
       this = ConstructEmptyMatrix_lsc(rows, columns)
    ENDIF

  END SUBROUTINE ConstructEmptyMatrixSub_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Internal only. Create a sparse matrix with a certain number of columns
  !! and rows. Will allocate storage for the outer values, nothing else.
  !! @param[out] this the matrix being created. It will have the outer
  !! index allocated, but nothing else.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE FUNCTION ConstructEmptyMatrix_lsc(rows, columns, zero_in) RESULT(this)
    !! Parameters
    TYPE(Matrix_lsc) :: this
    INTEGER, INTENT(IN) :: columns, rows
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    INCLUDE "sparse_includes/ConstructEmptyMatrix.f90"
  END FUNCTION ConstructEmptyMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ConstructMatrixFromFileSub_lsc(this, file_name)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name

    this = ConstructMatrixFromFile_lsc(file_name)
  END SUBROUTINE ConstructMatrixFromFileSub_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] this the matrix being constructed.
  !! @param[in] file_name name of the file.
  FUNCTION ConstructMatrixFromFile_lsc(file_name) RESULT(this)
    !! Parameters
    TYPE(Matrix_lsc) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Triplet_c) :: temporary
    REAL(NTREAL) :: real_val, comp_val

#define ISCOMPLEX
#include "sparse_includes/ConstructMatrixFromFile.f90"
#undef ISCOMPLEX

  END FUNCTION ConstructMatrixFromFile_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE ConstructMatrixFromTripletListSub_lsc(this, triplet_list, &
       & rows, columns)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(INOUT) :: this
    TYPE(TripletList_c), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN) :: rows, columns

    this = ConstructMatrixFromTripletList_lsc(triplet_list, rows, columns)

  END SUBROUTINE ConstructMatrixFromTripletListSub_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! The triplet list must be sorted to efficiently fill in the matrix. This
  !! constructor assumes \b you have already sorted the triplet list.
  !! @param[out] this the matrix being constructed
  !! @param[in] triplet_list a list of triplet values. They must be sorted.
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE FUNCTION ConstructMatrixFromTripletList_lsc(triplet_list,rows,columns) &
       & RESULT(this)
    !! Parameters
    TYPE(Matrix_lsc) :: this
    TYPE(TripletList_c), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN) :: rows, columns

    INCLUDE "sparse_includes/ConstructMatrixFromTripletList.f90"
  END FUNCTION ConstructMatrixFromTripletList_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  !! @param[inout] this the matrix to free up
  PURE SUBROUTINE DestructMatrix_lsc(this)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(INOUT) :: this

    INCLUDE "sparse_includes/DestructMatrix.f90"
  END SUBROUTINE DestructMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE CopyMatrix_lsc(matA, matB)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: matA
    TYPE(Matrix_lsc), INTENT(INOUT) :: matB

    INCLUDE "sparse_includes/CopyMatrix.f90"
  END SUBROUTINE CopyMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  !! @param[in] this the matrix.
  !! @result number of rows.
  PURE FUNCTION GetMatrixRows_lsc(this) RESULT(rows)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this

    INCLUDE "sparse_includes/GetMatrixRows.f90"
  END FUNCTION GetMatrixRows_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  !! @param[in] this the matrix.
  !! @result number of columns.
  PURE FUNCTION GetMatrixColumns_lsc(this) RESULT(columns)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this

    INCLUDE "sparse_includes/GetMatrixColumns.f90"
  END FUNCTION GetMatrixColumns_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixRow_lsc(this, row_number, row_out)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    TYPE(Matrix_lsc), INTENT(INOUT) :: row_out
    !! Temporary Variables
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: value_buffer

    INCLUDE "sparse_includes/ExtractMatrixRow.f90"
  END SUBROUTINE ExtractMatrixRow_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixColumn_lsc(this, column_number, column_out)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    TYPE(Matrix_lsc), INTENT(INOUT) :: column_out

    INCLUDE "sparse_includes/ExtractMatrixColumn.f90"
  END SUBROUTINE ExtractMatrixColumn_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! The current implementation has you go from matrix to triplet list,
  !! triplet list to transposed triplet list. The triplet list must then be
  !! sorted and then the return matrix is constructed.
  !! @param[in] this the matrix to be transposed.
  !! @param[out] matT the input matrix transposed.
  PURE SUBROUTINE TransposeMatrix_lsc(this, matT)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN)  :: this
    TYPE(Matrix_lsc), INTENT(INOUT) :: matT

    INCLUDE "sparse_includes/TransposeMatrix.f90"
  END SUBROUTINE TransposeMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !! to another.
  !! @param[in] mat_array 2d array of matrices to compose.
  !! @param[in] block_rows the number of rows of the array of blocks.
  !! @param[in] block_columns the number of columns of the array of blocks.
  !! @param[out] out_matrix the composed matrix.
  PURE SUBROUTINE ComposeMatrix_lsc(mat_array, block_rows, block_columns, &
       & out_matrix)
    !! Parameters
    TYPE(Matrix_lsc), DIMENSION(block_rows,block_columns), INTENT(IN) :: &
         & mat_array
    INTEGER, INTENT(IN) :: block_rows, block_columns
    TYPE(Matrix_lsc), INTENT(INOUT) :: out_matrix
    !! Local Data
    TYPE(Matrix_lsc), DIMENSION(block_columns) :: merged_columns
    TYPE(Matrix_lsc) :: Temp
    TYPE(Matrix_lsc), DIMENSION(block_rows,block_columns) :: mat_t

    INCLUDE "sparse_includes/ComposeMatrix.f90"
  END SUBROUTINE ComposeMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
  !! the first matrix are followed by the columns of the matrices in the list.
  !! @param[in] mat_list list of matrices to compose.
  !! @param[out] out_matrix = [Matrix 1 | Matrix 2, ...] .
  PURE SUBROUTINE ComposeMatrixColumns_lsc(mat_list, out_matrix)
    !! Parameters
    TYPE(Matrix_lsc), DIMENSION(:), INTENT(IN) :: mat_list
    TYPE(Matrix_lsc), INTENT(INOUT) :: out_matrix

    INCLUDE "sparse_includes/ComposeMatrixColumns.f90"
  END SUBROUTINE ComposeMatrixColumns_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  !! @param[in] this the matrix to split.
  !! @param[in] block_rows number of rows to split the matrix into.
  !! @param[in] block_columns number of columns to split the matrix into.
  !! @param[out] split_array a COLUMNxROW array for the output to go into.
  !! @param[in] block_size_row_in specifies the size of the  rows.
  !! @param[in] block_size_column_in specifies the size of the columns.
  SUBROUTINE SplitMatrix_lsc(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: block_rows, block_columns
    TYPE(Matrix_lsc), DIMENSION(:,:), INTENT(INOUT) :: split_array
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
    !! Local Data
    TYPE(Matrix_lsc), DIMENSION(block_columns) :: column_split
    TYPE(Matrix_lsc), DIMENSION(block_rows) :: row_split
    TYPE(Matrix_lsc) :: Temp

    INCLUDE "sparse_includes/SplitMatrix.f90"
  END SUBROUTINE SplitMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a matrix into into small blocks based on the specified offsets.
  !! @param[in] this matrix to perform this operation on.
  !! @param[in] num_blocks number of blocks to split into
  !! @param[out] block_sizes the sizes used for splitting.
  !! @param[out] split_list 1D array of blocks.
  PURE SUBROUTINE SplitMatrixColumns_lsc(this, num_blocks, block_sizes, &
       & split_list)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: num_blocks
    INTEGER, DIMENSION(num_blocks), INTENT(IN) :: block_sizes
    TYPE(Matrix_lsc), DIMENSION(num_blocks), INTENT(INOUT) :: split_list

    INCLUDE "sparse_includes/SplitMatrixColumns.f90"
  END SUBROUTINE SplitMatrixColumns_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list from a matrix.
  !! @param[in] this the matrix to construct the triplet list from.
  !! @param[out] triplet_list the triplet list we created.
  PURE SUBROUTINE MatrixToTripletList_lsc(this, triplet_list)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    !! Local Variables
    TYPE(Triplet_c) :: temporary

    INCLUDE "sparse_includes/MatrixToTripletList.f90"
  END SUBROUTINE MatrixToTripletList_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintMatrix_lsc(this, file_name_in)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TripletList_c) :: triplet_list

#define ISCOMPLEX
#include "sparse_includes/PrintMatrix.f90"
#undef ISCOMPLEX
  END SUBROUTINE PrintMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Every value in the matrix is changed into its complex conjugate.
  !! @param[inout] this the matrix to compute the complex conjugate of.
  PURE SUBROUTINE ConjugateMatrix_lsc(this)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(INOUT) :: this

    this%values = CONJG(this%values)
  END SUBROUTINE ConjugateMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a complex matrix to a real matrix.
  !! @param[in] cin the starting matrix.
  !! @param[out] rout real valued matrix.
  SUBROUTINE ConvertMatrixType_lsrtolsc(cin, rout)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN)    :: cin
    TYPE(Matrix_lsr), INTENT(INOUT) :: rout
    !! Local Variables
    TYPE(TripletList_c) :: in_list
    TYPE(TripletList_r) :: out_list

    CALL MatrixToTripletList(cin, in_list)
    CALL ConvertTripletListType(in_list, out_list)
    CALL ConstructMatrixFromTripletList(rout, out_list, cin%rows, cin%columns)

    CALL DestructTripletList(in_list)
    CALL DestructTripletList(out_list)

  END SUBROUTINE ConvertMatrixType_lsrtolsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a real matrix to a complex matrix.
  !! @param[in] rin the starting matrix.
  !! @param[out] cout complex valued matrix.
  SUBROUTINE ConvertMatrixType_lsctolsr(rin, cout)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)    :: rin
    TYPE(Matrix_lsc), INTENT(INOUT) :: cout
    !! Local Variables
    TYPE(TripletList_r) :: in_list
    TYPE(TripletList_c) :: out_list

    CALL MatrixToTripletList(rin, in_list)
    CALL ConvertTripletListType(in_list, out_list)
    CALL ConstructMatrixFromTripletList(cout, out_list, rin%rows, rin%columns)

    CALL DestructTripletList(in_list)
    CALL DestructTripletList(out_list)

  END SUBROUTINE ConvertMatrixType_lsctolsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixSModule
