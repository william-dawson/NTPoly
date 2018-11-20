!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Helper Routines for Computing The Cholesky Decomposition
MODULE CholeskyModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, MPINTINTEGER
  USE DMatrixModule, ONLY : Matrix_ldr
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceAndComposeMatrixSizes, &
       & ReduceAndComposeMatrixData, ReduceAndComposeMatrixCleanup, &
       & TestReduceSizeRequest, TestReduceInnerRequest, TestReduceDataRequest
  USE PSMatrixModule, ONLY : Matrix_ps, FillMatrixFromTripletList
  USE ProcessGridModule, ONLY : ProcessGrid_t
  USE SolverParametersModule, ONLY : SolverParameters_t
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, TransposeMatrix, &
       & DestructMatrix
  USE SVectorModule, ONLY : DotSparseVectors
  USE TripletListModule, ONLY : TripletList_r, AppendToTripletList, &
       & DestructTripletList
  USE TripletModule, ONLY : Triplet_r
  USE NTMPIModule
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
  PUBLIC :: CholeskyFactor_r
  PUBLIC :: CholeskyFactor_c
  PUBLIC :: CholeskyInvert_r
  PUBLIC :: CholeskyInvert_c
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
  !> Compute the cholesky factor of a real matrix (this = LLT).
  SUBROUTINE CholeskyFactor_r(mat, factor, parameters)
    !> The matrix to invert.
    TYPE(Matrix_ps), INTENT(IN)  :: mat
    !> The inverse of that matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: factor
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN) :: parameters
  END SUBROUTINE CholeskyFactor_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the cholesky factor of a complex matrix (this = LLT).
  SUBROUTINE CholeskyFactor_c(mat, factor, parameters)
    !> The matrix to invert.
    TYPE(Matrix_ps), INTENT(IN)  :: mat
    !> The inverse of that matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: factor
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN) :: parameters
  END SUBROUTINE CholeskyFactor_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Invert a lower triangular matrix (real).
  SUBROUTINE CholeskyInvert_r(mat, factor, parameters)
    !> The matrix to invert.
    TYPE(Matrix_ps), INTENT(IN)  :: mat
    !> The inverse of that matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: factor
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN) :: parameters
  END SUBROUTINE CholeskyInvert_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Invert a lower triangular matrix (complex).
  SUBROUTINE CholeskyInvert_c(mat, factor, parameters)
    !> The matrix to invert.
    TYPE(Matrix_ps), INTENT(IN)  :: mat
    !> The inverse of that matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: factor
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN) :: parameters
  END SUBROUTINE CholeskyInvert_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to insert a value into a sparse vector.
  PURE SUBROUTINE AppendToVector_r(values_per, indices, values, insert_row, &
       & insert_value)
    !> Values per row.
    INTEGER, INTENT(INOUT) :: values_per
    !> Indices associated with each value.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    !> Values.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    !> Row to insert into.
    INTEGER, INTENT(IN) :: insert_row
    !> Value to insert.
    REAL(NTREAL), INTENT(IN) :: insert_value

    INCLUDE "solver_includes/AppendToVector.f90"
  END SUBROUTINE AppendToVector_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to broadcast a sparse vector
  SUBROUTINE BroadcastVector_r(num_values, indices, values, root, comm)
    !> Number of values we are broadcasting.
    INTEGER, INTENT(INOUT) :: num_values
    !> Indices to broadcast.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    !> Values to broadcast.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    !> Root from which we broadcast.
    INTEGER, INTENT(IN) :: root
    !> Communicator to broadcast along.
    INTEGER, INTENT(INOUT) :: comm
    !! Local
    INTEGER :: err

    INCLUDE "solver_includes/BroadcastVector.f90"
    CALL MPI_Bcast(values(:num_values), num_values, MPINTREAL, root, comm, err)

  END SUBROUTINE BroadcastVector_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the vector holding the accumulated diagonal values
  SUBROUTINE ConstructDiag_r(AMat, process_grid, dense_a, diag)
    !> AMat the matrix we're working on (for meta data).
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> The process grid we're operating on.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !> A dense representation of the values.
    TYPE(Matrix_ldr), INTENT(IN) :: dense_a
    !> Diagonal values computed.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: diag

    INCLUDE "solver_includes/ConstructDiag.f90"
    CALL MPI_Allgatherv(MPI_IN_PLACE, diags_per_proc(process_grid%my_row+1), &
         & MPINTREAL, diag, diags_per_proc, diag_displ, MPINTREAL, &
         & process_grid%column_comm, ierr)
  END SUBROUTINE ConstructDiag_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a lookup for columns
  SUBROUTINE ConstructRankLookup(AMat, process_grid, col_root_lookup)
    !> Matrix we're computing.
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> Grid we're computing along.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !> The lookup we're computing.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: col_root_lookup
    !! Local Variables
    INTEGER, DIMENSION(process_grid%num_process_columns) :: cols_per_proc
    INTEGER, DIMENSION(process_grid%num_process_columns) :: d_cols_per_proc
    INTEGER :: II
    INTEGER :: ierr

    cols_per_proc(process_grid%my_column+1) = AMat%local_columns
    CALL MPI_Allgather(MPI_IN_PLACE, 1, MPINTINTEGER, cols_per_proc, 1, &
         & MPINTINTEGER, process_grid%row_comm, ierr)
    d_cols_per_proc(1) = 0
    DO II = 2, process_grid%num_process_columns
       d_cols_per_proc(II) = d_cols_per_proc(II-1) + cols_per_proc(II-1)
    END DO
    col_root_lookup(AMat%start_column:AMat%end_column - 1) = process_grid%row_rank
    CALL MPI_Allgatherv(MPI_IN_PLACE, AMat%local_columns, MPINTINTEGER, &
         & col_root_lookup, cols_per_proc, d_cols_per_proc, MPINTINTEGER, &
         & process_grid%row_comm, ierr)

  END SUBROUTINE ConstructRankLookup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine which computes sparse dot products across processors.
  !> Computes the dot product of one vector with several others.
  SUBROUTINE DotAllHelper_r(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, out_values, comm)
    !> The length of vector i.
    INTEGER, INTENT(IN) :: num_values_i
    !> Tn array with the length of vectors j.
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    !> The index value of the sparse vector i.
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    !> The indices of the vectors j.
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    !> The values of the sparse vector i.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    !> The values of the vectors j.
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    !> The dot product values for each vector j.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    !> The communicator to reduce along.
    INTEGER, INTENT(INOUT) :: comm

    INCLUDE "solver_includes/DotAllHelper.f90"

  END SUBROUTINE DotAllHelper_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine which computes sparse dot products across processors.
  !> Computes the dot product of one vector with several others.
  !> The pivoted version has the number of local pivots to work on as a
  !> parameter.
  SUBROUTINE DotAllPivoted_r(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, pivot_vector, num_local_pivots, out_values, comm)
    !> The length of vector i.
    INTEGER, INTENT(IN) :: num_values_i
    !> Tn array with the length of vectors j.
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    !> The index value of the sparse vector i.
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    !> The indices of the vectors j.
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    !> The values of the sparse vector i.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    !> The values of the vectors j.
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    !> Vector storing the pivot values.
    INTEGER, DIMENSION(:), INTENT(IN) :: pivot_vector
    !> Number of pivots.
    INTEGER, INTENT(IN) :: num_local_pivots
    !> The dot product values for each vector j.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    !> The communicator to reduce along.
    INTEGER, INTENT(INOUT) :: comm

    INCLUDE "solver_includes/DotAllPivoted.f90"

  END SUBROUTINE DotAllPivoted_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine that gathers the matrices in the same column into one.
  SUBROUTINE GatherMatrixColumn_r(local_matrix, column_matrix, process_grid)
    !> The local matrix on each process.
    TYPE(Matrix_lsr), INTENT(IN) :: local_matrix
    !> The final result.
    TYPE(Matrix_lsr), INTENT(INOUT) :: column_matrix
    !> The process grid to operate on.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !! Local Variables
    TYPE(Matrix_lsr) :: local_matrixT

    INCLUDE "solver_includes/GatherMatrixColumn.f90"
  END SUBROUTINE GatherMatrixColumn_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the pivot vector.
  SUBROUTINE GetPivot_r(AMat, process_grid, start_index, pivot_vector, diag, &
       & index, value, local_pivots, num_local_pivots)
    !> The matrix we're working on.
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> The process grid to compute on.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !> The current pivot vector.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: pivot_vector
    !> The diagonal values.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: diag
    !> The start index to look
    INTEGER, INTENT(IN) :: start_index
    !> The pivot index selected.
    INTEGER, INTENT(OUT) :: index
    !> The pivot value.
    REAL(NTREAL), INTENT(OUT) :: value
    !> The local pivot values to modify.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: local_pivots
    !> Number of pivots stored locally.
    INTEGER, INTENT(OUT) :: num_local_pivots
    !! Local Variables
    REAL(NTREAL) :: temp_diag
    DOUBLE PRECISION, DIMENSION(2) :: max_diag

    INCLUDE "solver_includes/GetPivot.f90"

  END SUBROUTINE GetPivot_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Unpack to a global matrix.
  SUBROUTINE UnpackCholesky_r(values_per_column, index, values, LMat)
    !> The number of values in a column.
    INTEGER, DIMENSION(:), INTENT(IN) :: values_per_column
    !> Index values.
    INTEGER, DIMENSION(:,:), INTENT(IN) :: index
    !> Actual values.
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values
    !> Matrix to unpack into.
    TYPE(Matrix_ps), INTENT(INOUT) :: LMat

    INCLUDE "solver_includes/UnpackCholesky.f90"
  END SUBROUTINE UnpackCholesky_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE CholeskyModule
