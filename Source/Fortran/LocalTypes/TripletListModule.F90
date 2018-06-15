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
  END TYPE TripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  TYPE, PUBLIC :: TripletList_c
     !> Internal representation of the data.
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: DATA
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
  END TYPE TripletList_c
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DATATYPE REAL(NTREAL)
#define MPIDATATYPE MPINTREAL
#define TLISTTYPE TripletList_r
#define TTYPE Triplet_r

#define ConstructTripletList ConstructTripletList_r
#define DestructTripletList DestructTripletList_r
#define ResizeTripletList ResizeTripletList_r
#define AppendToTripletList AppendToTripletList_r
#define AccumulateTripletList AccumulateTripletList_r
#define SetTripletAt SetTripletAt_r
#define GetTripletAt GetTripletAt_r
#define SortTripletList SortTripletList_r
#define SymmetrizeTripletList SymmetrizeTripletList_r
#define GetTripletListSize GetTripletListSize_r
#define RedistributeTripletLists RedistributeTripletLists_r
#define ShiftTripletList ShiftTripletList_r
#define SortDenseTripletList SortDenseTripletList_r

#include "includes/TripletListImpl.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] triplet_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList(triplet_list, pattern_type)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT)  :: triplet_list
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(TTYPE) :: temporary, temporary_transpose
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
  END SUBROUTINE SymmetrizeTripletList

#undef ConstructTripletList
#undef DestructTripletList
#undef ResizeTripletList
#undef AppendToTripletList
#undef AccumulateTripletList
#undef SetTripletAt
#undef GetTripletAt
#undef SortTripletList
#undef SymmetrizeTripletList
#undef GetTripletListSize
#undef RedistributeTripletLists
#undef ShiftTripletList
#undef SortDenseTripletList

#undef TTYPE
#undef TLISTTYPE
#undef MPIDATATYPE
#undef DATATYPE

#define DATATYPE COMPLEX(NTCOMPLEX)
#define MPIDATATYPE MPINTCOMPLEX
#define TLISTTYPE TripletList_c
#define TTYPE Triplet_c

#define ConstructTripletList ConstructTripletList_c
#define DestructTripletList DestructTripletList_c
#define ResizeTripletList ResizeTripletList_c
#define AppendToTripletList AppendToTripletList_c
#define AccumulateTripletList AccumulateTripletList_c
#define SetTripletAt SetTripletAt_c
#define GetTripletAt GetTripletAt_c
#define SortTripletList SortTripletList_c
#define SymmetrizeTripletList SymmetrizeTripletList_c
#define GetTripletListSize GetTripletListSize_c
#define RedistributeTripletLists RedistributeTripletLists_c
#define ShiftTripletList ShiftTripletList_c
#define SortDenseTripletList SortDenseTripletList_c

#include "includes/TripletListImpl.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] triplet_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList(triplet_list, pattern_type)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT)  :: triplet_list
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(TTYPE) :: temporary, temporary_transpose
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
  END SUBROUTINE SymmetrizeTripletList

#undef ConstructTripletList
#undef DestructTripletList
#undef ResizeTripletList
#undef AppendToTripletList
#undef AccumulateTripletList
#undef SetTripletAt
#undef GetTripletAt
#undef SortTripletList
#undef SymmetrizeTripletList
#undef GetTripletListSize
#undef RedistributeTripletLists
#undef ShiftTripletList
#undef SortDenseTripletList

#undef TTYPE
#undef TLISTTYPE
#undef MPIDATATYPE
#undef DATATYPE

END MODULE TripletListModule
