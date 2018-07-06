!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Helper Routines for Computing The Cholesky Decomposition
MODULE CholeskyModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DMatrixModule, ONLY : Matrix_ldr
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceMatrixSizes, &
       & ReduceAndComposeMatrixData, ReduceAndComposeMatrixCleanup
  USE PSMatrixModule, ONLY : Matrix_ps, FillMatrixFromTripletList
  USE ProcessGridModule, ONLY : ProcessGrid_t
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, TransposeMatrix, &
       & DestructMatrix
  USE SVectorModule, ONLY : DotSparseVectors
  USE TripletListModule, ONLY : TripletList_r, AppendToTripletList, &
       & DestructTripletList
  USE TripletModule, ONLY : Triplet_r
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: AppendToVector
  PUBLIC :: BroadcastVector
  PUBLIC :: ConstructDiag
  PUBLIC :: ConstructRankLookup
  PUBLIC :: DotAllHelper
  PUBLIC :: DotAllPivoted
  PUBLIC :: GatherMatrixColumn
  PUBLIC :: GetPivot
  PUBLIC :: UnpackCholesky
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE AppendToVector
     MODULE PROCEDURE AppendToVector_r
  END INTERFACE
  INTERFACE BroadcastVector
     MODULE PROCEDURE BroadcastVector_r
  END INTERFACE
  INTERFACE ConstructDiag
     MODULE PROCEDURE ConstructDiag_r
  END INTERFACE
  INTERFACE DotAllHelper
     MODULE PROCEDURE DotAllHelper_r
  END INTERFACE
  INTERFACE DotAllPivoted
     MODULE PROCEDURE DotAllPivoted_r
  END INTERFACE
  INTERFACE GatherMatrixColumn
     MODULE PROCEDURE GatherMatrixColumn_r
  END INTERFACE
  INTERFACE GetPivot
     MODULE PROCEDURE GetPivot_r
  END INTERFACE
  INTERFACE UnpackCholesky
     MODULE PROCEDURE UnpackCholesky_r
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to insert a value into a sparse vector.
  PURE SUBROUTINE AppendToVector_r(values_per, indices, values, insert_row, &
       & insert_value)
    !! Parameters
    INTEGER, INTENT(INOUT) :: values_per
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    INTEGER, INTENT(IN) :: insert_row
    REAL(NTREAL), INTENT(IN) :: insert_value

    INCLUDE "includes/AppendToVector.F90"
  END SUBROUTINE AppendToVector_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to broadcast a sparse vector
  SUBROUTINE BroadcastVector_r(num_values, indices, values, root, comm)
    !! Parameters
    INTEGER, INTENT(INOUT) :: num_values
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    INTEGER, INTENT(IN) :: root
    INTEGER, INTENT(INOUT) :: comm
    !! Local
    INTEGER :: err

    INCLUDE "includes/BroadcastVector.F90"
    CALL MPI_Bcast(values(:num_values), num_values, MPINTREAL, root, comm, err)

  END SUBROUTINE BroadcastVector_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the vector holding the accumulated diagonal values
  SUBROUTINE ConstructDiag_r(AMat, process_grid, dense_a, diag)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    TYPE(Matrix_ldr), INTENT(IN) :: dense_a
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: diag

    INCLUDE "includes/ConstructDiag.F90"
    CALL MPI_Allgatherv(MPI_IN_PLACE, diags_per_proc(process_grid%my_row+1), &
         & MPINTREAL, diag, diags_per_proc, diag_displ, MPINTREAL, &
         & process_grid%column_comm, ierr)
  END SUBROUTINE ConstructDiag_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a lookup for columns
  SUBROUTINE ConstructRankLookup(AMat, process_grid, col_root_lookup)
    !! Root Lookups
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    INTEGER, DIMENSION(:), INTENT(INOUT) :: col_root_lookup
    !! Local Variables
    INTEGER, DIMENSION(process_grid%num_process_columns) :: cols_per_proc
    INTEGER, DIMENSION(process_grid%num_process_columns) :: d_cols_per_proc
    INTEGER :: II
    INTEGER :: ierr

    cols_per_proc(process_grid%my_column+1) = AMat%local_columns
    CALL MPI_Allgather(MPI_IN_PLACE, 1, MPI_INTEGER, cols_per_proc, 1, &
         & MPI_INTEGER, process_grid%row_comm, ierr)
    d_cols_per_proc(1) = 0
    DO II = 2, process_grid%num_process_columns
       d_cols_per_proc(II) = d_cols_per_proc(II-1) + cols_per_proc(II-1)
    END DO
    col_root_lookup(AMat%start_column:AMat%end_column - 1) = process_grid%row_rank
    CALL MPI_Allgatherv(MPI_IN_PLACE, AMat%local_columns, MPI_INTEGER, &
         & col_root_lookup, cols_per_proc, d_cols_per_proc, MPI_INTEGER, &
         & process_grid%row_comm, ierr)

  END SUBROUTINE ConstructRankLookup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine which computes sparse dot products across processors.
  !! Computes the dot product of one vector with several others.
  !! @param[in] num_values_i the length of vector i.
  !! @param[in] indices_i the index value of the sparse vector i.
  !! @param[in] values_i the values of the sparse vector i.
  !! @param[in] num_values_j an array with the length of vectors j.
  !! @param[in] indices_j the indices of the vectors j.
  !! @param[in] values_j the values of the vectors j.
  !! @param[out] out_values the dot product values for each vector j.
  !! @param[in] comm the communicator to reduce along.
  SUBROUTINE DotAllHelper_r(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, out_values, comm)
    !! Parameters
    INTEGER, INTENT(IN) :: num_values_i
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    INTEGER, INTENT(INOUT) :: comm

    INCLUDE "includes/DotAllHelper.F90"

  END SUBROUTINE DotAllHelper_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DotAllPivoted_r(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, pivot_vector, num_local_pivots, out_values, comm)
    !! Parameters
    INTEGER, INTENT(IN) :: num_values_i
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    INTEGER, DIMENSION(:), INTENT(IN) :: pivot_vector
    INTEGER, INTENT(IN) :: num_local_pivots
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    INTEGER, INTENT(INOUT) :: comm

    INCLUDE "includes/DotAllPivoted.F90"

  END SUBROUTINE DotAllPivoted_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine that gathers the matrices in the same column into one.
  !! @param[in] local_matrix the local matrix on each process.
  !! @param[out] column_matrix the final result.
  SUBROUTINE GatherMatrixColumn_r(local_matrix, column_matrix, process_grid)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: local_matrix
    TYPE(Matrix_lsr), INTENT(INOUT) :: column_matrix
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !! Local Variables
    TYPE(Matrix_lsr) :: local_matrixT

    INCLUDE "includes/GatherMatrixColumn.F90"
  END SUBROUTINE GatherMatrixColumn_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetPivot_r(AMat, process_grid, start_index, pivot_vector, diag, &
       & index, value, local_pivots, num_local_pivots)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    INTEGER, DIMENSION(:), INTENT(INOUT) :: pivot_vector
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: diag
    INTEGER, INTENT(IN) :: start_index
    INTEGER, INTENT(OUT) :: index
    REAL(NTREAL), INTENT(OUT) :: value
    INTEGER, DIMENSION(:), INTENT(INOUT) :: local_pivots
    INTEGER, INTENT(OUT) :: num_local_pivots
    !! Local Variables
    REAL(NTREAL) :: temp_diag
    DOUBLE PRECISION, DIMENSION(2) :: max_diag

    INCLUDE "includes/GetPivot.F90"

  END SUBROUTINE GetPivot_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE UnpackCholesky_r(values_per_column, index, values, LMat)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN) :: values_per_column
    INTEGER, DIMENSION(:,:), INTENT(IN) :: index
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values
    TYPE(Matrix_ps), INTENT(INOUT) :: LMat

    INCLUDE "includes/UnpackCholesky.F90"
  END SUBROUTINE UnpackCholesky_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE CholeskyModule
