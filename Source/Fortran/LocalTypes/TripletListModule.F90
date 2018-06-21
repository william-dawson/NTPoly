!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Lists of Triplets.
MODULE TripletListModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE TripletModule, ONLY : Triplet, Triplet_r, Triplet_c
  USE MatrixMarketModule, ONLY : MM_SYMMETRIC, MM_SKEW_SYMMETRIC, MM_HERMITIAN
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE ISO_C_BINDING, ONLY : c_int
  USE MPI
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SortTripletList
  PUBLIC :: RedistributeTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, ABSTRACT, PUBLIC :: TripletList
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
   CONTAINS
     !! Basic Routines
     PROCEDURE(Construct_base), DEFERRED :: Init
     PROCEDURE(Destruct_base), DEFERRED :: Destruct
     PROCEDURE(Resize_base), DEFERRED :: Resize
     PROCEDURE :: GetSize => GetTripletListSize_b
     !! Manipulate Elements
     PROCEDURE(Append_base), DEFERRED :: Append
     PROCEDURE(Set_base), DEFERRED :: Set
     PROCEDURE(Get_base), DEFERRED :: Get
     !! Global Operations
     PROCEDURE(Symmetrize_base), DEFERRED :: Symmetrize
     PROCEDURE(Shift_base), DEFERRED :: Shift
  END TYPE TripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(TripletList), PUBLIC :: TripletList_r
     !> Internal representation of the data.
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: data
   CONTAINS
     !! Basic Routines
     PROCEDURE :: Init => ConstructTripletList_r
     PROCEDURE :: Destruct => DestructTripletList_r
     PROCEDURE :: Resize => ResizeTripletList_r
     !! Manipulate Elements
     PROCEDURE :: Append => AppendToTripletList_r
     PROCEDURE :: Set => SetTripletAt_r
     PROCEDURE :: Get => GetTripletAt_r
     !! Global Operations
     PROCEDURE :: Symmetrize => SymmetrizeTripletList_r
     PROCEDURE :: Shift => ShiftTripletList_r
  END TYPE TripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(TripletList), PUBLIC :: TripletList_c
     !> Internal representation of the data.
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: data
   CONTAINS
     !! Basic Routines
     PROCEDURE :: Init => ConstructTripletList_c
     PROCEDURE :: Destruct => DestructTripletList_c
     PROCEDURE :: Resize => ResizeTripletList_c
     !! Manipulate Elements
     PROCEDURE :: Append => AppendToTripletList_c
     PROCEDURE :: Set => SetTripletAt_c
     PROCEDURE :: Get => GetTripletAt_c
     !! Global Operations
     PROCEDURE :: Symmetrize => SymmetrizeTripletList_c
     PROCEDURE :: Shift => ShiftTripletList_c
  END TYPE TripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE SortTripletList
    MODULE PROCEDURE SortTripletList_r
    MODULE PROCEDURE SortTripletList_c
  END INTERFACE
  INTERFACE RedistributeTripletList
    MODULE PROCEDURE RedistributeTripletLists_r
    MODULE PROCEDURE RedistributeTripletLists_c
  END INTERFACE
  INTERFACE SortDenseTripletList
    MODULE PROCEDURE SortDenseTripletList_r
    MODULE PROCEDURE SortDenseTripletList_c
  END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ABSTRACT INTERFACE
    SUBROUTINE Construct_base(this, size_in)
      USE ISO_C_BINDING, ONLY : c_int
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(INOUT) :: this
      INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in
    END SUBROUTINE Construct_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Destruct_base(this)
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(INOUT) :: this
    END SUBROUTINE Destruct_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Resize_base(this, size)
      USE ISO_C_BINDING, ONLY : c_int
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(INOUT) :: this
      INTEGER(KIND=c_int), INTENT(IN) :: size
    END SUBROUTINE Resize_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Append_base(this, triplet_value)
      USE TripletModule, ONLY : Triplet
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(INOUT) :: this
      CLASS(Triplet), INTENT(IN)        :: triplet_value
    END SUBROUTINE Append_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Set_base(this, index, triplet_value)
      USE ISO_C_BINDING, ONLY : c_int
      USE TripletModule, ONLY : Triplet
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(INOUT) :: this
      INTEGER(kind=c_int), INTENT(IN)  :: index
      CLASS(Triplet), INTENT(IN)        :: triplet_value
    END SUBROUTINE Set_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Get_base(this, index, triplet_value)
      USE ISO_C_BINDING, ONLY : c_int
      USE TripletModule, ONLY : Triplet
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(IN)  :: this
      INTEGER(kind=c_int), INTENT(IN) :: index
      CLASS(Triplet), INTENT(INOUT)   :: triplet_value
    END SUBROUTINE Get_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Symmetrize_base(this, pattern_type)
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)               :: pattern_type
    END SUBROUTINE Symmetrize_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Shift_base(this, row_shift, column_shift)
      IMPORT :: TripletList
      IMPLICIT NONE
      CLASS(TripletList), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: row_shift
      INTEGER, INTENT(IN) :: column_shift
    END SUBROUTINE Shift_base
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  !! @param[inout] this the triplet list to construct.
  !! @param[in] size_in the length of the triplet list (optional, default=0).
  PURE SUBROUTINE ConstructTripletList_r(this, size_in)
    !! Parameters
    CLASS(TripletList_r), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    INCLUDE "triplet_includes/ConstructTripletList.f90"

  END SUBROUTINE ConstructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  !! @param[inout] this the triplet list to construct.
  !! @param[in] size_in the length of the triplet list (optional, default=0).
  PURE SUBROUTINE ConstructTripletList_c(this, size_in)
    !! Parameters
    CLASS(TripletList_c), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    INCLUDE "triplet_includes/ConstructTripletList.f90"

  END SUBROUTINE ConstructTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] this the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList_r(this)
    !! Parameters
    CLASS(TripletList_r), INTENT(INOUT) :: this

    INCLUDE "triplet_includes/DestructTripletList.f90"

  END SUBROUTINE DestructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] this the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList_c(this)
    !! Parameters
    CLASS(TripletList_c), INTENT(INOUT) :: this

    INCLUDE "triplet_includes/DestructTripletList.f90"

  END SUBROUTINE DestructTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  !! @param[inout] this the triplet list to resize.
  !! @param[in] size to resize to.
  PURE SUBROUTINE ResizeTripletList_r(this, size)
    !! Parameters
    CLASS(TripletList_r), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: temporary_data

    INCLUDE "triplet_includes/ResizeTripletList.f90"

  END SUBROUTINE ResizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  !! @param[inout] this the triplet list to resize.
  !! @param[in] size to resize to.
  PURE SUBROUTINE ResizeTripletList_c(this, size)
    !! Parameters
    CLASS(TripletList_c), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: temporary_data

    INCLUDE "triplet_includes/ResizeTripletList.f90"

  END SUBROUTINE ResizeTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  !! @param[in] triplet_list list to get the size of.
  !! @return list_size the number of entries in the triplet list.
  PURE FUNCTION GetTripletListSize_b(triplet_list) RESULT(list_size)
    !! Parameters
    CLASS(TripletList), INTENT(IN)  :: triplet_list
    INTEGER :: list_size

    list_size = triplet_list%CurrentSize

  END FUNCTION GetTripletListSize_b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  !! @param[inout] this the triplet list to append to.
  !! @param[in] triplet_value the value to append.
  PURE SUBROUTINE AppendToTripletList_r(this, triplet_value)
    !! Parameters
    CLASS(TripletList_r), INTENT(INOUT) :: this
    CLASS(Triplet), INTENT(IN)          :: triplet_value
    !! Local data
    INTEGER :: new_size
    INTEGER :: old_size

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_r)
      INCLUDE "triplet_includes/AppendToTripletList.f90"
    END SELECT

  END SUBROUTINE AppendToTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  !! @param[inout] this the triplet list to append to.
  !! @param[in] triplet_value the value to append.
  PURE SUBROUTINE AppendToTripletList_c(this, triplet_value)
    !! Parameters
    CLASS(TripletList_c), INTENT(INOUT) :: this
    CLASS(Triplet), INTENT(IN)        :: triplet_value
    !! Local data
    INTEGER :: new_size
    INTEGER :: old_size

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_c)
      INCLUDE "triplet_includes/AppendToTripletList.f90"
    END SELECT

  END SUBROUTINE AppendToTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  !! @param[inout] this the triplet list to set.
  !! @param[in] index the index at which to set the triplet.
  !! @param[in] triplet_value the value of the triplet to set.
  PURE SUBROUTINE SetTripletAt_r(this,index,triplet_value)
    !! Parameters
    CLASS(TripletList_r), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN)     :: index
    CLASS(Triplet), INTENT(IN)        :: triplet_value

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_r)
      this%data(index) = triplet_value
    END SELECT

  END SUBROUTINE SetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  !! @param[inout] this the triplet list to set.
  !! @param[in] index the index at which to set the triplet.
  !! @param[in] triplet_value the value of the triplet to set.
  PURE SUBROUTINE SetTripletAt_c(this,index,triplet_value)
    !! Parameters
    CLASS(TripletList_c), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN)     :: index
    CLASS(Triplet), INTENT(IN)        :: triplet_value

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_c)
      this%data(index) = triplet_value
    END SELECT

  END SUBROUTINE SetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  !! @param[in] this the triplet list to get the value from.
  !! @param[in] index the index from which to get the triplet.
  !! @param[out] triplet_value the extracted triplet value.
  PURE SUBROUTINE GetTripletAt_r(this,index,triplet_value)
    !! Parameters
    CLASS(TripletList_r), INTENT(IN) :: this
    INTEGER(kind=c_int), INTENT(IN)  :: index
    CLASS(Triplet), INTENT(INOUT)  :: triplet_value

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_r)
      triplet_value%index_row = this%data(index)%index_row
      triplet_value%index_column = this%data(index)%index_column
      triplet_value%point_value = this%data(index)%point_value
    END SELECT

  END SUBROUTINE GetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  !! @param[in] this the triplet list to get the value from.
  !! @param[in] index the index from which to get the triplet.
  !! @param[out] triplet_value the extracted triplet value.
  PURE SUBROUTINE GetTripletAt_c(this,index,triplet_value)
    !! Parameters
    CLASS(TripletList_c), INTENT(IN) :: this
    INTEGER(kind=c_int), INTENT(IN)  :: index
    CLASS(Triplet), INTENT(INOUT)  :: triplet_value

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_c)
      triplet_value%index_row = this%data(index)%index_row
      triplet_value%index_column = this%data(index)%index_column
      triplet_value%point_value = this%data(index)%point_value
    END SELECT

  END SUBROUTINE GetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] triplet_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList_r(this, pattern_type)
    !! Parameters
    CLASS(TripletList_r), INTENT(INOUT)  :: this
    INTEGER, INTENT(IN) :: pattern_type
    !! Local Variables
    TYPE(Triplet_r) :: temporary

    INCLUDE "triplet_includes/SymmetrizeTripletList.f90"
  END SUBROUTINE SymmetrizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] triplet_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList_c(this, pattern_type)
    !! Parameters
    CLASS(TripletList_c), INTENT(INOUT)  :: this
    INTEGER, INTENT(IN) :: pattern_type
    !! Local Variables
    TYPE(Triplet_c) :: temporary

    INCLUDE "triplet_includes/SymmetrizeTripletList.f90"
  END SUBROUTINE SymmetrizeTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !! Frequently, we have a triplet list that comes from the global matrix which
  !! we would like to shift into a local matrix. In that case, just pass
  !! the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_r(this, row_shift, column_shift)
    !! Parameters
    CLASS(TripletList_r), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: row_shift
    INTEGER, INTENT(IN) :: column_shift

    INCLUDE "triplet_includes/ShiftTripletList.f90"

  END SUBROUTINE ShiftTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !! Frequently, we have a triplet list that comes from the global matrix which
  !! we would like to shift into a local matrix. In that case, just pass
  !! the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_c(this, row_shift, column_shift)
    !! Parameters
    CLASS(TripletList_c), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: row_shift
    INTEGER, INTENT(IN) :: column_shift

    INCLUDE "triplet_includes/ShiftTripletList.f90"

  END SUBROUTINE ShiftTripletList_c
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
    CLASS(TripletList_r), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList_r), INTENT(OUT) :: sorted_list
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_r) :: temporary

    INCLUDE "triplet_includes/SortTripletList.f90"

  END SUBROUTINE SortTripletList_r
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
    CLASS(TripletList_c), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList_c), INTENT(OUT) :: sorted_list
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_c) :: temporary

    INCLUDE "triplet_includes/SortTripletList.f90"

  END SUBROUTINE SortTripletList_c
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

    INCLUDE "triplet_includes/SortDenseTripletList.f90"

  END SUBROUTINE SortDenseTripletList_r
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

    INCLUDE "triplet_includes/SortDenseTripletList.f90"

  END SUBROUTINE SortDenseTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !! Takes in a list of triplet lists, one list for each processor. Then the
  !! all to all redistribution is performed along the given communicator.
  !! @param[in] triplet_lists a list of triplet lists, one for each process.
  !! @param[inout] comm the mpi communicator to redistribute along.
  !! @param[out] local_data_out the resulting local triplet list.
  SUBROUTINE RedistributeTripletLists_r(triplet_lists, comm, local_data_out)
    !! Parameters
    CLASS(TripletList_r), DIMENSION(:), INTENT(IN) :: triplet_lists
    INTEGER, INTENT(INOUT) :: comm
    CLASS(TripletList_r), INTENT(INOUT) :: local_data_out
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
  !! Takes in a list of triplet lists, one list for each processor. Then the
  !! all to all redistribution is performed along the given communicator.
  !! @param[in] triplet_lists a list of triplet lists, one for each process.
  !! @param[inout] comm the mpi communicator to redistribute along.
  !! @param[out] local_data_out the resulting local triplet list.
  SUBROUTINE RedistributeTripletLists_c(triplet_lists, comm, local_data_out)
    !! Parameters
    CLASS(TripletList_c), DIMENSION(:), INTENT(IN) :: triplet_lists
    INTEGER, INTENT(INOUT) :: comm
    CLASS(TripletList_c), INTENT(INOUT) :: local_data_out
    !! Local data (type specific)
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_c) :: temp_triplet

#define MPIDATATYPE MPINTCOMPLEX
#include "triplet_includes/RedistributeTripletLists.f90"
#undef MPIDATATYPE

  END SUBROUTINE RedistributeTripletLists_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletListModule
