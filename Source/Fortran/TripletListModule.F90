!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Lists of Triplets.
MODULE TripletListModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX, &
       & MPINTINTEGER
  USE TripletModule, ONLY : Triplet_r, Triplet_c, CompareTriplets, &
       & ConvertTripletType
  USE MatrixMarketModule, ONLY : MM_SYMMETRIC, MM_SKEW_SYMMETRIC, MM_HERMITIAN
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  TYPE :: TripletList_r
     !> Internal representation of the data.
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: DATA
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
  END TYPE TripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  TYPE :: TripletList_c
     !> Internal representation of the data.
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: DATA
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
  END TYPE TripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TripletList_r
  PUBLIC :: TripletList_c
  PUBLIC :: ConstructTripletList
  PUBLIC :: CopyTripletList
  PUBLIC :: DestructTripletList
  PUBLIC :: ResizeTripletList
  PUBLIC :: AppendToTripletList
  PUBLIC :: SetTripletAt
  PUBLIC :: GetTripletAt
  PUBLIC :: SortTripletList
  PUBLIC :: SymmetrizeTripletList
  PUBLIC :: GetTripletListSize
  PUBLIC :: RedistributeTripletLists
  PUBLIC :: AllGatherTripletList
  PUBLIC :: ShiftTripletList
  PUBLIC :: ConvertTripletListType
  PUBLIC :: MergeTripletLists
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructTripletList
     MODULE PROCEDURE ConstructTripletListSup_r
     MODULE PROCEDURE ConstructTripletListSup_c
  END INTERFACE ConstructTripletList
  INTERFACE CopyTripletList
     MODULE PROCEDURE CopyTripletList_r
     MODULE PROCEDURE CopyTripletList_c
     MODULE PROCEDURE CopyTripletList_rc
  END INTERFACE CopyTripletList
  INTERFACE DestructTripletList
     MODULE PROCEDURE DestructTripletList_r
     MODULE PROCEDURE DestructTripletList_c
  END INTERFACE DestructTripletList
  INTERFACE ResizeTripletList
     MODULE PROCEDURE ResizeTripletList_r
     MODULE PROCEDURE ResizeTripletList_c
  END INTERFACE ResizeTripletList
  INTERFACE AppendToTripletList
     MODULE PROCEDURE AppendToTripletList_r
     MODULE PROCEDURE AppendToTripletList_c
  END INTERFACE AppendToTripletList
  INTERFACE SetTripletAt
     MODULE PROCEDURE SetTripletAt_r
     MODULE PROCEDURE SetTripletAt_c
  END INTERFACE SetTripletAt
  INTERFACE GetTripletAt
     MODULE PROCEDURE GetTripletAt_r
     MODULE PROCEDURE GetTripletAt_c
  END INTERFACE GetTripletAt
  INTERFACE SortTripletList
     MODULE PROCEDURE SortTripletList_r
     MODULE PROCEDURE SortTripletList_c
  END INTERFACE SortTripletList
  INTERFACE SortDenseTripletList
     MODULE PROCEDURE SortDenseTripletList_r
     MODULE PROCEDURE SortDenseTripletList_c
  END INTERFACE SortDenseTripletList
  INTERFACE SymmetrizeTripletList
     MODULE PROCEDURE SymmetrizeTripletList_r
     MODULE PROCEDURE SymmetrizeTripletList_c
  END INTERFACE SymmetrizeTripletList
  INTERFACE GetTripletListSize
     MODULE PROCEDURE GetTripletListSize_r
     MODULE PROCEDURE GetTripletListSize_c
  END INTERFACE GetTripletListSize
  INTERFACE RedistributeTripletLists
     MODULE PROCEDURE RedistributeTripletLists_r
     MODULE PROCEDURE RedistributeTripletLists_c
  END INTERFACE RedistributeTripletLists
  INTERFACE AllGatherTripletList
     MODULE PROCEDURE AllGatherTripletList_r
     MODULE PROCEDURE AllGatherTripletList_c
  END INTERFACE AllGatherTripletList
  INTERFACE ShiftTripletList
     MODULE PROCEDURE ShiftTripletList_r
     MODULE PROCEDURE ShiftTripletList_c
  END INTERFACE ShiftTripletList
  INTERFACE ConvertTripletListType
     MODULE PROCEDURE ConvertTripletListToReal
     MODULE PROCEDURE ConvertTripletListToComplex
  END INTERFACE ConvertTripletListType
  INTERFACE MergeTripletLists
     MODULE PROCEDURE MergeTripletLists_r
     MODULE PROCEDURE MergeTripletLists_c
  END INTERFACE MergeTripletLists
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for constructing a triplet list.
  PURE SUBROUTINE ConstructTripletListSup_r(this, size_in)
    !> The triplet list to construct.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> The length of the triplet list (default = 0).
    INTEGER, INTENT(IN), OPTIONAL :: size_in

#include "triplet_includes/ConstructTripletList.f90"

  END SUBROUTINE ConstructTripletListSup_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for constructing a triplet list.
  PURE SUBROUTINE ConstructTripletListSup_c(this, size_in)
    !> The triplet list to construct.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> The length of the triplet list (default = 0).
    INTEGER, INTENT(IN), OPTIONAL :: size_in

#include "triplet_includes/ConstructTripletList.f90"

  END SUBROUTINE ConstructTripletListSup_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  PURE SUBROUTINE DestructTripletList_r(this)
    !> The triplet list to destruct.
    TYPE(TripletList_r), INTENT(INOUT) :: this

#include "triplet_includes/DestructTripletList.f90"

  END SUBROUTINE DestructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  PURE SUBROUTINE DestructTripletList_c(this)
    !> The triplet list to destruct.
    TYPE(TripletList_c), INTENT(INOUT) :: this

#include "triplet_includes/DestructTripletList.f90"

  END SUBROUTINE DestructTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a triplet list (real).
  SUBROUTINE CopyTripletList_r(tripA, tripB)
    !> The triplet list to copy.
    TYPE(TripletList_r), INTENT(IN) :: tripA
    !> tripB = tripA
    TYPE(TripletList_r), INTENT(INOUT) :: tripB

#include "triplet_includes/CopyTripletList.f90"
  END SUBROUTINE CopyTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a triplet list (complex).
  SUBROUTINE CopyTripletList_c(tripA, tripB)
    !> The triplet list to copy.
    TYPE(TripletList_c), INTENT(IN) :: tripA
    !> tripB = tripA
    TYPE(TripletList_c), INTENT(INOUT) :: tripB

#include "triplet_includes/CopyTripletList.f90"
  END SUBROUTINE CopyTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy and upcast a triplet list (real -> complex).
  SUBROUTINE CopyTripletList_rc(tripA, tripB)
    !> The triplet list to copy.
    TYPE(TripletList_r), INTENT(IN) :: tripA
    !> tripB = tripA
    TYPE(TripletList_c), INTENT(INOUT) :: tripB
    !! Local varaibles
    INTEGER II

    CALL ConstructTripletList(tripB, tripA%CurrentSize)
    DO II = 1, tripA%CurrentSize
       tripB%DATA(II)%index_row = tripA%DATA(II)%index_row
       tripB%DATA(II)%index_column = tripA%DATA(II)%index_column
       tripB%DATA(II)%point_value = &
           & CMPLX(tripA%DATA(II)%point_value, KIND=NTCOMPLEX)
    END DO
  END SUBROUTINE CopyTripletList_rc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  PURE SUBROUTINE ResizeTripletList_r(this, size)
    !> The triplet list to resize.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> Size to resize to.
    INTEGER, INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: temporary_data

#include "triplet_includes/ResizeTripletList.f90"

  END SUBROUTINE ResizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  PURE SUBROUTINE ResizeTripletList_c(this, size)
    !> The triplet list to resize.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> Size to resize to.
    INTEGER, INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: temporary_data

#include "triplet_includes/ResizeTripletList.f90"

  END SUBROUTINE ResizeTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  PURE SUBROUTINE AppendToTripletList_r(this, triplet_value)
    !> This the triplet list to append to.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> The value to append.
    TYPE(Triplet_r), INTENT(IN)        :: triplet_value

#include "triplet_includes/AppendToTripletList.f90"

  END SUBROUTINE AppendToTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  PURE SUBROUTINE AppendToTripletList_c(this, triplet_value)
    !> This the triplet list to append to.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> The value to append.
    TYPE(Triplet_c), INTENT(IN)        :: triplet_value

#include "triplet_includes/AppendToTripletList.f90"

  END SUBROUTINE AppendToTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  PURE SUBROUTINE SetTripletAt_r(this,index,triplet_value)
    !> The triplet list to set.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> The index at which to set the triplet.
    INTEGER, INTENT(IN)    :: index
    !> The value of the triplet to set.
    TYPE(Triplet_r), INTENT(IN)        :: triplet_value

#include "triplet_includes/SetTripletAt.f90"
  END SUBROUTINE SetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  PURE SUBROUTINE SetTripletAt_c(this,index,triplet_value)
    !> The triplet list to set.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> The index at which to set the triplet.
    INTEGER, INTENT(IN)    :: index
    !> The value of the triplet to set.
    TYPE(Triplet_c), INTENT(IN)        :: triplet_value

#include "triplet_includes/SetTripletAt.f90"
  END SUBROUTINE SetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  PURE SUBROUTINE GetTripletAt_r(this,index,triplet_value)
    !> The triplet list to get the value from.
    TYPE(TripletList_r), INTENT(IN) :: this
    !> The index from which to get the triplet.
    INTEGER, INTENT(IN) :: index
    !> The extracted triplet value.
    TYPE(Triplet_r), INTENT(OUT)    :: triplet_value

#include "triplet_includes/GetTripletAt.f90"
  END SUBROUTINE GetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  PURE SUBROUTINE GetTripletAt_c(this,index,triplet_value)
    !> The triplet list to get the value from.
    TYPE(TripletList_c), INTENT(IN) :: this
    !> The index from which to get the triplet.
    INTEGER, INTENT(IN) :: index
    !> The extracted triplet value.
    TYPE(Triplet_c), INTENT(OUT)    :: triplet_value

#include "triplet_includes/GetTripletAt.f90"
  END SUBROUTINE GetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !> Implementation is based on bucket sort. This is why it needs the number of
  !> matrix columns. Bubble sort is used within a bucket.
  PURE SUBROUTINE SortTripletList_r(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !> List to be sorted.
    TYPE(TripletList_r), INTENT(IN)  :: input_list
    !> This is the highest column value in the list.
    INTEGER, INTENT(IN) :: matrix_columns
    !> This is the highest row value in the list.
    INTEGER, INTENT(IN) :: matrix_rows
    !> A now sorted version of the list. This routine will allocate it.
    TYPE(TripletList_r), INTENT(OUT) :: sorted_list
    !> False if you do not need the final bubble sort.
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_r) :: trip

#include "triplet_includes/SortTripletList.f90"

  END SUBROUTINE SortTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !> Implementation is based on bucket sort. This is why it needs the number of
  !> matrix columns. Bubble sort is used within a bucket.
  PURE SUBROUTINE SortTripletList_c(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !> List to be sorted.
    TYPE(TripletList_c), INTENT(IN)  :: input_list
    !> This is the highest column value in the list.
    INTEGER, INTENT(IN) :: matrix_columns
    !> This is the highest row value in the list.
    INTEGER, INTENT(IN) :: matrix_rows
    !> A now sorted version of the list. This routine will allocate it.
    TYPE(TripletList_c), INTENT(OUT) :: sorted_list
    !> False if you do not need the final bubble sort.
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_c) :: trip

#include "triplet_includes/SortTripletList.f90"

  END SUBROUTINE SortTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  PURE FUNCTION GetTripletListSize_r(triplet_list) RESULT(list_size)
    !> List to get the size of.
    TYPE(TripletList_r), INTENT(IN)  :: triplet_list
    !> The number of entries in the triplet list.
    INTEGER :: list_size

#include "triplet_includes/GetTripletListSize.f90"

  END FUNCTION GetTripletListSize_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  PURE FUNCTION GetTripletListSize_c(triplet_list) RESULT(list_size)
    !> List to get the size of.
    TYPE(TripletList_c), INTENT(IN)  :: triplet_list
    !> The number of entries in the triplet list.
    INTEGER :: list_size

#include "triplet_includes/GetTripletListSize.f90"

  END FUNCTION GetTripletListSize_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !> Takes in a list of triplet lists, one list for each processor. Then the
  !> all to all redistribution is performed along the given communicator.
  SUBROUTINE RedistributeTripletLists_r(triplet_lists, comm, local_data_out)
    !> A list of triplet lists, one for each process.
    TYPE(TripletList_r), DIMENSION(:), INTENT(IN) :: triplet_lists
    !> The mpi communicator to redistribute along.
    INTEGER, INTENT(IN) :: comm
    !> The resulting local triplet list.
    TYPE(TripletList_r), INTENT(INOUT) :: local_data_out
    !! Local data (type specific)
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_r) :: temp_triplet

#define MPIDATATYPE MPINTREAL
#include "triplet_includes/RedistributeTripletLists.f90"
#undef MPIDATATYPE

  END SUBROUTINE RedistributeTripletLists_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !> Takes in a list of triplet lists, one list for each processor. Then the
  !> all to all redistribution is performed along the given communicator.
  SUBROUTINE RedistributeTripletLists_c(triplet_lists, comm, local_data_out)
    !> A list of triplet lists, one for each process.
    TYPE(TripletList_c), DIMENSION(:), INTENT(IN) :: triplet_lists
    !> The mpi communicator to redistribute along.
    INTEGER, INTENT(IN) :: comm
    !> The resulting local triplet list.
    TYPE(TripletList_c), INTENT(INOUT) :: local_data_out
    !! Local data (type specific)
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_c) :: temp_triplet

#define MPIDATATYPE MPINTCOMPLEX
#include "triplet_includes/RedistributeTripletLists.f90"
#undef MPIDATATYPE

  END SUBROUTINE RedistributeTripletLists_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gather triplet lists from a set of processors.
  SUBROUTINE AllGatherTripletList_r(triplet_in, comm, gathered_out)
    !> Locally held triplet list
    TYPE(TripletList_r), INTENT(IN) :: triplet_in
    !> The mpi communicator to gather along.
    INTEGER, INTENT(IN) :: comm
    !> The resulting gathered triplet list.
    TYPE(TripletList_r), INTENT(INOUT) :: gathered_out
    !! Local data (type specific)
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_r) :: temp_triplet

#define MPIDATATYPE MPINTREAL
#include "triplet_includes/GatherTripletList.f90"
#undef MPIDATATYPE

  END SUBROUTINE AllGatherTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gather triplet lists from a set of processors.
  SUBROUTINE AllGatherTripletList_c(triplet_in, comm, gathered_out)
    !> Locally held triplet list
    TYPE(TripletList_c), INTENT(IN) :: triplet_in
    !> The mpi communicator to gather along.
    INTEGER, INTENT(IN) :: comm
    !> The resulting gathered triplet list.
    TYPE(TripletList_c), INTENT(INOUT) :: gathered_out
    !! Local data (type specific)
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_c) :: temp_triplet

#define MPIDATATYPE MPINTCOMPLEX
#include "triplet_includes/GatherTripletList.f90"
#undef MPIDATATYPE

  END SUBROUTINE AllGatherTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !> Frequently, we have a triplet list that comes from the global matrix which
  !> we would like to shift into a local matrix. In that case, just pass
  !> the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_r(triplet_list, row_shift, column_shift)
    !> The triplet list to shift.
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    !> The row offset to shift by.
    INTEGER, INTENT(IN) :: row_shift
    !> The column offset to shift by.
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: II

#include "triplet_includes/ShiftTripletList.f90"

  END SUBROUTINE ShiftTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !> Frequently, we have a triplet list that comes from the global matrix which
  !> we would like to shift into a local matrix. In that case, just pass
  !> the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_c(triplet_list, row_shift, column_shift)
    !> The triplet list to shift.
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    !> The row offset to shift by.
    INTEGER, INTENT(IN) :: row_shift
    !> The column offset to shift by.
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: II

#include "triplet_includes/ShiftTripletList.f90"

  END SUBROUTINE ShiftTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !> dense.
  PURE SUBROUTINE SortDenseTripletList_r(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !> The list to sort.
    TYPE(TripletList_r), INTENT(IN)  :: input_list
    !> Number of columns for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_columns
    !> Number of rows for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_rows
    !> Sorted and ready to use for building matrices.
    TYPE(TripletList_r), INTENT(OUT) :: sorted_list
    !! Local Variables
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: value_buffer

#include "triplet_includes/SortDenseTripletList.f90"

  END SUBROUTINE SortDenseTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !> dense.
  PURE SUBROUTINE SortDenseTripletList_c(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !> The list to sort.
    TYPE(TripletList_c), INTENT(IN)  :: input_list
    !> Number of columns for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_columns
    !> Number of rows for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_rows
    !> Sorted and ready to use for building matrices.
    TYPE(TripletList_c), INTENT(OUT) :: sorted_list
    !! Local Variables
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: value_buffer

#include "triplet_includes/SortDenseTripletList.f90"

  END SUBROUTINE SortDenseTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !> symmetry type.
  SUBROUTINE SymmetrizeTripletList_r(triplet_list, pattern_type)
    !> List to be symmetrized.
    TYPE(TripletList_r), INTENT(INOUT)  :: triplet_list
    !> Type of symmetry.
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(Triplet_r) :: trip, trip_t
    INTEGER :: II
    INTEGER :: initial_size

    initial_size = triplet_list%CurrentSize
    SELECT CASE(pattern_type)
    CASE(MM_SYMMETRIC)
       DO II = 1, initial_size
          CALL GetTripletAt(triplet_list, II, trip)
          IF (trip%index_column .NE. trip%index_row) THEN
             trip_t%index_row = trip%index_column
             trip_t%index_column = trip%index_row
             trip_t%point_value = trip%point_value
             CALL AppendToTripletList(triplet_list, trip_t)
          END IF
       END DO
    CASE(MM_SKEW_SYMMETRIC)
       DO II = 1, initial_size
          CALL GetTripletAt(triplet_list, II, trip)
          IF (trip%index_column .NE. trip%index_row) THEN
             trip_t%index_row = trip%index_column
             trip_t%index_column = trip%index_row
             trip_t%point_value = -1.0 * trip%point_value
             CALL AppendToTripletList(triplet_list, trip_t)
          END IF
       END DO
    END SELECT
  END SUBROUTINE SymmetrizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !> symmetry type.
  SUBROUTINE SymmetrizeTripletList_c(triplet_list, pattern_type)
    !> List to be symmetrized.
    TYPE(TripletList_c), INTENT(INOUT)  :: triplet_list
    !> Type of symmetry.
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(Triplet_c) :: trip, trip_t
    INTEGER :: II
    INTEGER :: initial_size

    initial_size = triplet_list%CurrentSize
    SELECT CASE(pattern_type)
    CASE(MM_SYMMETRIC)
       DO II = 1, initial_size
          CALL GetTripletAt(triplet_list, II, trip)
          IF (trip%index_column .NE. trip%index_row) THEN
             trip_t%index_row = trip%index_column
             trip_t%index_column = trip%index_row
             trip_t%point_value = trip%point_value
             CALL AppendToTripletList(triplet_list, trip_t)
          END IF
       END DO
    CASE(MM_HERMITIAN)
       DO II = 1, initial_size
          CALL GetTripletAt(triplet_list, II, trip)
          IF (trip%index_column .NE. trip%index_row) THEN
             trip_t%index_row = trip%index_column
             trip_t%index_column = trip%index_row
             trip_t%point_value = CONJG(trip%point_value)
             CALL AppendToTripletList(triplet_list, trip_t)
          END IF
       END DO
    CASE(MM_SKEW_SYMMETRIC)
       DO II = 1, initial_size
          CALL GetTripletAt(triplet_list, II, trip)
          IF (trip%index_column .NE. trip%index_row) THEN
             trip_t%index_row = trip%index_column
             trip_t%index_column = trip%index_row
             trip_t%point_value = -1.0*trip%point_value
             CALL AppendToTripletList(triplet_list, trip_t)
          END IF
       END DO
    END SELECT
  END SUBROUTINE SymmetrizeTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a complex triplet list to a real triplet list.
  SUBROUTINE ConvertTripletListToReal(cin_triplet, rout_triplet)
    !> The starting triplet list.
    TYPE(TripletList_c), INTENT(IN)    :: cin_triplet
    !> Real valued triplet list.
    TYPE(TripletList_r), INTENT(INOUT) :: rout_triplet
    !! Local Variables
    INTEGER :: II

    CALL ConstructTripletList(rout_triplet, cin_triplet%CurrentSize)
    DO II = 1, cin_triplet%CurrentSize
       CALL ConvertTripletType(cin_triplet%DATA(II), rout_triplet%DATA(II))
    END DO

  END SUBROUTINE ConvertTripletListToReal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a real triplet to a complex triplet list.
  SUBROUTINE ConvertTripletListToComplex(rin_triplet, cout_triplet)
    !> The starting triplet list.
    TYPE(TripletList_r), INTENT(IN)    :: rin_triplet
    !> Complex valued triplet list.
    TYPE(TripletList_c), INTENT(INOUT) :: cout_triplet
    !! Local Variables
    INTEGER :: II

    CALL ConstructTripletList(cout_triplet, rin_triplet%CurrentSize)
    DO II = 1, rin_triplet%CurrentSize
       CALL ConvertTripletType(rin_triplet%DATA(II), cout_triplet%DATA(II))
    END DO

  END SUBROUTINE ConvertTripletListToComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Merge a list of TripletList objects together.
  SUBROUTINE MergeTripletLists_r(this, tlists)
    !> The starting triplet list. This will be completely replaced.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> An array of triplet list objects.
    TYPE(TripletList_r), DIMENSION(:), INTENT(IN) :: tlists

#include "triplet_includes/MergeTripletLists.f90"

  END SUBROUTINE MergeTripletLists_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Merge a list of TripletList objects together.
  SUBROUTINE MergeTripletLists_c(this, tlists)
    !> The starting triplet list. This will be completely replaced.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> An array of triplet list objects.
    TYPE(TripletList_c), DIMENSION(:), INTENT(IN) :: tlists

#include "triplet_includes/MergeTripletLists.f90"

  END SUBROUTINE MergeTripletLists_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletListModule
