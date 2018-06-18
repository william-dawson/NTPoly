!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Lists of Triplets.
MODULE TripletListModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE TripletModule, ONLY : Triplet_r, Triplet_c, CompareTriplets
  USE MatrixMarketModule, ONLY : MM_SYMMETRIC, MM_SKEW_SYMMETRIC, MM_HERMITIAN
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE ISO_C_BINDING, ONLY : c_int
  USE MPI
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TripletList_r
  PUBLIC :: TripletList_c
  PUBLIC :: DestructTripletList
  PUBLIC :: ResizeTripletList
  PUBLIC :: AppendToTripletList
  PUBLIC :: AccumulateTripletList
  PUBLIC :: SetTripletAt
  PUBLIC :: GetTripletAt
  PUBLIC :: SortTripletList
  PUBLIC :: SymmetrizeTripletList
  PUBLIC :: GetTripletListSize
  PUBLIC :: RedistributeTripletLists
  PUBLIC :: ShiftTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE TripletList_r
     MODULE PROCEDURE ConstructTripletList_r
  END INTERFACE
  INTERFACE TripletList_c
     MODULE PROCEDURE ConstructTripletList_c
  END INTERFACE
  INTERFACE ConstructTripletList
     MODULE PROCEDURE ConstructTripletListSup_r
     MODULE PROCEDURE ConstructTripletListSup_c
  END INTERFACE
  INTERFACE DestructTripletList
     MODULE PROCEDURE DestructTripletList_r
     MODULE PROCEDURE DestructTripletList_c
  END INTERFACE
  INTERFACE ResizeTripletList
     MODULE PROCEDURE ResizeTripletList_r
     MODULE PROCEDURE ResizeTripletList_c
  END INTERFACE
  INTERFACE AppendToTripletList
     MODULE PROCEDURE AppendToTripletList_r
     MODULE PROCEDURE AppendToTripletList_c
  END INTERFACE
  INTERFACE AccumulateTripletList
     MODULE PROCEDURE AccumulateTripletList_r
     MODULE PROCEDURE AccumulateTripletList_c
  END INTERFACE
  INTERFACE SetTripletAt
     MODULE PROCEDURE SetTripletAt_r
     MODULE PROCEDURE SetTripletAt_c
  END INTERFACE
  INTERFACE GetTripletAt
     MODULE PROCEDURE GetTripletAt_r
     MODULE PROCEDURE GetTripletAt_c
  END INTERFACE
  INTERFACE SortTripletList
     MODULE PROCEDURE SortTripletList_r
     MODULE PROCEDURE SortTripletList_c
  END INTERFACE
  INTERFACE SortDenseTripletList
     MODULE PROCEDURE SortDenseTripletList_r
     MODULE PROCEDURE SortDenseTripletList_c
  END INTERFACE
  INTERFACE SymmetrizeTripletList
     MODULE PROCEDURE SymmetrizeTripletList_r
     MODULE PROCEDURE SymmetrizeTripletList_c
  END INTERFACE
  INTERFACE GetTripletListSize
     MODULE PROCEDURE GetTripletListSize_r
     MODULE PROCEDURE GetTripletListSize_c
  END INTERFACE
  INTERFACE RedistributeTripletLists
     MODULE PROCEDURE RedistributeTripletLists_r
     MODULE PROCEDURE RedistributeTripletLists_c
  END INTERFACE
  INTERFACE ShiftTripletList
     MODULE PROCEDURE ShiftTripletList_r
     MODULE PROCEDURE ShiftTripletList_c
  END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  TYPE, PUBLIC :: TripletList_r
     !> Internal representation of the data.
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: DATA
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
   CONTAINS
     FINAL :: DestructTripletList_r
  END TYPE TripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  TYPE, PUBLIC :: TripletList_c
     !> Internal representation of the data.
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: DATA
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
   CONTAINS
     FINAL :: DestructTripletList_c
  END TYPE TripletList_c
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE ConstructTripletListSup_r(this, size_in)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    IF (PRESENT(size_in)) THEN
       this = ConstructTripletList_r(size_in)
    ELSE
       this = ConstructTripletList_r()
    END IF
  END SUBROUTINE ConstructTripletListSup_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  !! @param[inout] this the triplet list to construct.
  !! @param[in] size_in the length of the triplet list (optional, default=0).
  PURE FUNCTION ConstructTripletList_r(size_in) RESULT(this)
    !! Parameters
    TYPE(TripletList_r) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    INCLUDE "triplet_includes/ConstructTripletList.f90"

  END FUNCTION ConstructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] this the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList_r(this)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT) :: this

    INCLUDE "triplet_includes/DestructTripletList.f90"

  END SUBROUTINE DestructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  !! @param[inout] this the triplet list to resize.
  !! @param[in] size to resize to.
  PURE SUBROUTINE ResizeTripletList_r(this, size)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: temporary_data

    INCLUDE "triplet_includes/ResizeTripletList.f90"

  END SUBROUTINE ResizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  !! @param[inout] this the triplet list to append to.
  !! @param[in] triplet_value the value to append.
  PURE SUBROUTINE AppendToTripletList_r(this, triplet_value)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT) :: this
    TYPE(Triplet_r), INTENT(IN)        :: triplet_value

    INCLUDE "triplet_includes/AppendToTripletList.f90"

  END SUBROUTINE AppendToTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> (Just for a related project)
  PURE SUBROUTINE AccumulateTripletList_r(this, triplet_value)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT) :: this
    TYPE(Triplet_r), INTENT(IN)        :: triplet_value

    INCLUDE "triplet_includes/AccumulateTripletList.f90"

  END SUBROUTINE AccumulateTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  !! @param[inout] this the triplet list to set.
  !! @param[in] index the index at which to set the triplet.
  !! @param[in] triplet_value the value of the triplet to set.
  PURE SUBROUTINE SetTripletAt_r(this,index,triplet_value)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN)    :: index
    TYPE(Triplet_r), INTENT(IN)        :: triplet_value

    INCLUDE "triplet_includes/SetTripletAt.f90"
  END SUBROUTINE SetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  !! @param[in] this the triplet list to get the value from.
  !! @param[in] index the index from which to get the triplet.
  !! @param[out] triplet_value the extracted triplet value.
  PURE SUBROUTINE GetTripletAt_r(this,index,triplet_value)
    !! Parameters
    TYPE(TripletList_r), INTENT(IN) :: this
    INTEGER(kind=c_int), INTENT(IN) :: index
    TYPE(Triplet_r), INTENT(OUT)    :: triplet_value

    INCLUDE "triplet_includes/GetTripletAt.f90"
  END SUBROUTINE GetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !! Implementation is based on bucket sort. This is why it needs the number of
  !! matrix columns. Bubble sort is used within a bucket.
  !! @param[in] input_list list to be sorted.
  !! @param[in] matrix_columns this is the highest column value in the list.
  !! @param[in] matrix_rows this is the highest row value in the list.
  !! @param[in] bubble_in false if you don't need the final bubble sort.
  !! @param[out] sorted_list a now sorted version of the list. This routine
  !! will allocate it.
  PURE SUBROUTINE SortTripletList_r(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !! Parameters
    TYPE(TripletList_r), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList_r), INTENT(OUT) :: sorted_list
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_r) :: temporary

#include "triplet_includes/SortTripletList.f90"

  END SUBROUTINE SortTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  !! @param[in] triplet_list list to get the size of.
  !! @return list_size the number of entries in the triplet list.
  PURE FUNCTION GetTripletListSize_r(triplet_list) RESULT(list_size)
    !! Parameters
    TYPE(TripletList_r), INTENT(IN)  :: triplet_list
    INTEGER :: list_size

    INCLUDE "triplet_includes/GetTripletListSize.f90"

  END FUNCTION GetTripletListSize_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !! Takes in a list of triplet lists, one list for each processor. Then the
  !! all to all redistribution is performed along the given communicator.
  !! @param[in] triplet_lists a list of triplet lists, one for each process.
  !! @param[inout] comm the mpi communicator to redistribute along.
  !! @param[out] local_data_out the resulting local triplet list.
  SUBROUTINE RedistributeTripletLists_r(triplet_lists, comm, local_data_out)
    !! Parameters
    TYPE(TripletList_r), DIMENSION(:), INTENT(IN) :: triplet_lists
    INTEGER, INTENT(INOUT) :: comm
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
  !> Shift the rows and columns of a triplet list by set values.
  !! Frequently, we have a triplet list that comes from the global matrix which
  !! we would like to shift into a local matrix. In that case, just pass
  !! the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_r(triplet_list, row_shift, column_shift)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    INTEGER, INTENT(IN) :: row_shift
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: counter

    INCLUDE "triplet_includes/ShiftTripletList.f90"

  END SUBROUTINE ShiftTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !! dense.
  !! @param[in] input_list the list
  !! @param[in] matrix_columns for the corresponding matrix.
  !! @param[in] matrix_rows for the corresponding matrix.
  !! @param[out] sorted_list sorted and ready to use for building matrices.
  PURE SUBROUTINE SortDenseTripletList_r(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !! Parameters
    TYPE(TripletList_r), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList_r), INTENT(OUT) :: sorted_list
    !! Local Variables
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: value_buffer

#include "triplet_includes/SortDenseTripletList.f90"

  END SUBROUTINE SortDenseTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] triplet_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList_r(triplet_list, pattern_type)
    !! Parameters
    TYPE(TripletList_r), INTENT(INOUT)  :: triplet_list
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(Triplet_r) :: temporary, temporary_transpose
    INTEGER :: counter
    INTEGER :: initial_size

    initial_size = triplet_list%CurrentSize
    SELECT CASE(pattern_type)
    CASE(MM_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    CASE(MM_SKEW_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = -1.0*temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    END SELECT
  END SUBROUTINE SymmetrizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE ConstructTripletListSup_c(this, size_in)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    IF (PRESENT(size_in)) THEN
       this = ConstructTripletList_c(size_in)
    ELSE
       this = ConstructTripletList_c()
    END IF
  END SUBROUTINE ConstructTripletListSup_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  !! @param[inout] this the triplet list to construct.
  !! @param[in] size_in the length of the triplet list (optional, default=0).
  PURE FUNCTION ConstructTripletList_c(size_in) RESULT(this)
    !! Parameters
    TYPE(TripletList_c) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    INCLUDE "triplet_includes/ConstructTripletList.f90"

  END FUNCTION ConstructTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] this the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList_c(this)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT) :: this

    INCLUDE "triplet_includes/DestructTripletList.f90"

  END SUBROUTINE DestructTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  !! @param[inout] this the triplet list to resize.
  !! @param[in] size to resize to.
  PURE SUBROUTINE ResizeTripletList_c(this, size)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: temporary_data

    INCLUDE "triplet_includes/ResizeTripletList.f90"

  END SUBROUTINE ResizeTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  !! @param[inout] this the triplet list to append to.
  !! @param[in] triplet_value the value to append.
  PURE SUBROUTINE AppendToTripletList_c(this, triplet_value)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT) :: this
    TYPE(Triplet_c), INTENT(IN)        :: triplet_value

    INCLUDE "triplet_includes/AppendToTripletList.f90"

  END SUBROUTINE AppendToTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> (Just for a related project)
  PURE SUBROUTINE AccumulateTripletList_c(this, triplet_value)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT) :: this
    TYPE(Triplet_c), INTENT(IN)        :: triplet_value

    INCLUDE "triplet_includes/AccumulateTripletList.f90"

  END SUBROUTINE AccumulateTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  !! @param[inout] this the triplet list to set.
  !! @param[in] index the index at which to set the triplet.
  !! @param[in] triplet_value the value of the triplet to set.
  PURE SUBROUTINE SetTripletAt_c(this,index,triplet_value)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN)    :: index
    TYPE(Triplet_c), INTENT(IN)        :: triplet_value

    INCLUDE "triplet_includes/SetTripletAt.f90"
  END SUBROUTINE SetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  !! @param[in] this the triplet list to get the value from.
  !! @param[in] index the index from which to get the triplet.
  !! @param[out] triplet_value the extracted triplet value.
  PURE SUBROUTINE GetTripletAt_c(this,index,triplet_value)
    !! Parameters
    TYPE(TripletList_c), INTENT(IN) :: this
    INTEGER(kind=c_int), INTENT(IN) :: index
    TYPE(Triplet_c), INTENT(OUT)    :: triplet_value

    INCLUDE "triplet_includes/GetTripletAt.f90"
  END SUBROUTINE GetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !! Implementation is based on bucket sort. This is why it needs the number of
  !! matrix columns. Bubble sort is used within a bucket.
  !! @param[in] input_list list to be sorted.
  !! @param[in] matrix_columns this is the highest column value in the list.
  !! @param[in] matrix_rows this is the highest row value in the list.
  !! @param[in] bubble_in false if you don't need the final bubble sort.
  !! @param[out] sorted_list a now sorted version of the list. This routine
  !! will allocate it.
  PURE SUBROUTINE SortTripletList_c(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !! Parameters
    TYPE(TripletList_c), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList_c), INTENT(OUT) :: sorted_list
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_c) :: temporary

#include "triplet_includes/SortTripletList.f90"

  END SUBROUTINE SortTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  !! @param[in] triplet_list list to get the size of.
  !! @return list_size the number of entries in the triplet list.
  PURE FUNCTION GetTripletListSize_c(triplet_list) RESULT(list_size)
    !! Parameters
    TYPE(TripletList_c), INTENT(IN)  :: triplet_list
    INTEGER :: list_size

    INCLUDE "triplet_includes/GetTripletListSize.f90"

  END FUNCTION GetTripletListSize_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !! Takes in a list of triplet lists, one list for each processor. Then the
  !! all to all redistribution is performed along the given communicator.
  !! @param[in] triplet_lists a list of triplet lists, one for each process.
  !! @param[inout] comm the mpi communicator to redistribute along.
  !! @param[out] local_data_out the resulting local triplet list.
  SUBROUTINE RedistributeTripletLists_c(triplet_lists, comm, local_data_out)
    !! Parameters
    TYPE(TripletList_c), DIMENSION(:), INTENT(IN) :: triplet_lists
    INTEGER, INTENT(INOUT) :: comm
    TYPE(TripletList_c), INTENT(INOUT) :: local_data_out
    !! Local data (type specific)
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_c) :: temp_triplet

#define MPIDATATYPE MPINTREAL
#include "triplet_includes/RedistributeTripletLists.f90"
#undef MPIDATATYPE

  END SUBROUTINE RedistributeTripletLists_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !! Frequently, we have a triplet list that comes from the global matrix which
  !! we would like to shift into a local matrix. In that case, just pass
  !! the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_c(triplet_list, row_shift, column_shift)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    INTEGER, INTENT(IN) :: row_shift
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: counter

    INCLUDE "triplet_includes/ShiftTripletList.f90"

  END SUBROUTINE ShiftTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !! dense.
  !! @param[in] input_list the list
  !! @param[in] matrix_columns for the corresponding matrix.
  !! @param[in] matrix_rows for the corresponding matrix.
  !! @param[out] sorted_list sorted and ready to use for building matrices.
  PURE SUBROUTINE SortDenseTripletList_c(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !! Parameters
    TYPE(TripletList_c), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList_c), INTENT(OUT) :: sorted_list
    !! Local Variables
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: value_buffer

#include "triplet_includes/SortDenseTripletList.f90"

  END SUBROUTINE SortDenseTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] triplet_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList_c(triplet_list, pattern_type)
    !! Parameters
    TYPE(TripletList_c), INTENT(INOUT)  :: triplet_list
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(Triplet_c) :: temporary, temporary_transpose
    INTEGER :: counter
    INTEGER :: initial_size

    initial_size = triplet_list%CurrentSize
    SELECT CASE(pattern_type)
    CASE(MM_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    CASE(MM_HERMITIAN)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = CONJG(temporary%point_value)
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    CASE(MM_SKEW_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = -1.0*temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    END SELECT
  END SUBROUTINE SymmetrizeTripletList_c
END MODULE TripletListModule
