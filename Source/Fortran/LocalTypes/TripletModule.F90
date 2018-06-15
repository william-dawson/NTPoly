!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Value.
MODULE TripletModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE ISO_C_BINDING, ONLY : c_int
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetTriplet
  PUBLIC :: GetTripletValues
  PUBLIC :: CompareTriplets
  PUBLIC :: GetMPITripletType_r
  PUBLIC :: GetMPITripletType_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE SetTriplet
     MODULE PROCEDURE SetTriplet_r
     MODULE PROCEDURE SetTriplet_c
  END INTERFACE
  INTERFACE GetTripletValues
     MODULE PROCEDURE GetTripletValues_r
     MODULE PROCEDURE GetTripletValues_c
  END INTERFACE
  INTERFACE CompareTriplets
     MODULE PROCEDURE CompareTriplets_r
     MODULE PROCEDURE CompareTriplets_c
  END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a triplet of integer, integer, double.
  !! As this is related to matrix multiplication, the referencing indices are
  !! rows and columns.
  TYPE, PUBLIC :: Triplet_r
     INTEGER(kind=c_int)    :: index_column !< column value.
     INTEGER(kind=c_int)    :: index_row    !< row value.
     REAL(NTREAL) :: point_value  !< actual value at those indices.
  END TYPE Triplet_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a triplet of integer, integer, complex.
  !! As this is related to matrix multiplication, the referencing indices are
  !! rows and columns.
  TYPE, PUBLIC :: Triplet_c
     INTEGER(kind=c_int)    :: index_column !< column value.
     INTEGER(kind=c_int)    :: index_row    !< row value.
     COMPLEX(NTCOMPLEX) :: point_value  !< actual value at those indices.
  END TYPE Triplet_c
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DATATYPE REAL(NTREAL)
#define MPIDATATYPE MPINTREAL

#define INNERTYPE Triplet_r
#define SetTriplet SetTriplet_r
#define GetTripletValues GetTripletValues_r
#define CompareTriplets CompareTriplets_r
#define GetMPITripletType GetMPITripletType_r

#include "includes/TripletImpl.f90"

#undef SetTriplet
#undef GetTripletValues
#undef CompareTriplets
#undef GetMPITripletType

#undef INNERTYPE
#undef MPIDATATYPE
#undef DATATYPE

#define DATATYPE COMPLEX(NTCOMPLEX)
#define MPIDATATYPE MPINTCOMPLEX

#define INNERTYPE Triplet_c
#define SetTriplet SetTriplet_c
#define GetTripletValues GetTripletValues_c
#define CompareTriplets CompareTriplets_c
#define GetMPITripletType GetMPITripletType_c

#include "includes/TripletImpl.f90"

#undef SetTriplet
#undef GetTripletValues
#undef CompareTriplets
#undef GetMPITripletType

#undef INNERTYPE
#undef MPIDATATYPE
#undef DATATYPE
#undef make_name

END MODULE TripletModule
