!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for the triplet data type. 
!> Each one stores two indices and a value. This is related to sparse matrices, 
!> the referencing indices are usually rows and columns.
MODULE TripletModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX, &
       & MPINTINTEGER
  USE ErrorModule, ONLY : Error_t, CheckMPIError
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a triplet of integer, integer, double.
  TYPE, PUBLIC :: Triplet_r
     INTEGER :: index_column !< column value.
     INTEGER :: index_row    !< row value.
     REAL(NTREAL) :: point_value  !< actual value at those indices.
  END TYPE Triplet_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a triplet of integer, integer, complex.
  TYPE, PUBLIC :: Triplet_c
     INTEGER :: index_column !< column value.
     INTEGER :: index_row    !< row value.
     COMPLEX(NTCOMPLEX) :: point_value  !< actual value at those indices.
  END TYPE Triplet_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetTriplet
  PUBLIC :: GetTripletValues
  PUBLIC :: CompareTriplets
  PUBLIC :: GetMPITripletType_r
  PUBLIC :: GetMPITripletType_c
  PUBLIC :: ConvertTripletType
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE SetTriplet
     MODULE PROCEDURE SetTriplet_r
     MODULE PROCEDURE SetTriplet_c
  END INTERFACE SetTriplet
  INTERFACE GetTripletValues
     MODULE PROCEDURE GetTripletValues_r
     MODULE PROCEDURE GetTripletValues_c
  END INTERFACE GetTripletValues
  INTERFACE CompareTriplets
     MODULE PROCEDURE CompareTriplets_r
     MODULE PROCEDURE CompareTriplets_c
  END INTERFACE CompareTriplets
  INTERFACE ConvertTripletType
     MODULE PROCEDURE ConvertTripletToReal
     MODULE PROCEDURE ConvertTripletToComplex
  END INTERFACE ConvertTripletType
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet (real).
  PURE SUBROUTINE SetTriplet_r(this, index_column, index_row, point_value)
    !> The triplet to set the values of.
    TYPE(Triplet_r), INTENT(INOUT) :: this
    !> The column value.
    INTEGER, INTENT(IN) :: index_column
    !> The row value.
    INTEGER, INTENT(IN) :: index_row
    !> The value at that point.
    REAL(NTREAL), INTENT(IN) :: point_value

#include "triplet_includes/SetTriplet.f90"

  END SUBROUTINE SetTriplet_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet (complex).
  PURE SUBROUTINE SetTriplet_c(this, index_column, index_row, point_value)
    !> The triplet to set the values of.
    TYPE(Triplet_c), INTENT(INOUT) :: this
    !> The column value.
    INTEGER, INTENT(IN) :: index_column
    !> The row value.
    INTEGER, INTENT(IN) :: index_row
    !> The value at that point.
    COMPLEX(NTCOMPLEX), INTENT(IN) :: point_value

#include "triplet_includes/SetTriplet.f90"

  END SUBROUTINE SetTriplet_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  PURE SUBROUTINE GetTripletValues_r(this, index_column, index_row, &
       & point_value)
    !> The triplet to extract the values of.
    TYPE(Triplet_r), INTENT(IN) :: this
    !> Column value.
    INTEGER, INTENT(OUT) :: index_column
    !> Row value.
    INTEGER, INTENT(OUT) :: index_row
    !> Actual stored value.
    REAL(NTREAL), INTENT(OUT) :: point_value

#include "triplet_includes/GetTriplet.f90"

  END SUBROUTINE GetTripletValues_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  PURE SUBROUTINE GetTripletValues_c(this, index_column, index_row, &
       & point_value)
    !> The triplet to extract the values of.
    TYPE(Triplet_c), INTENT(IN) :: this
    !> Column value.
    INTEGER, INTENT(OUT) :: index_column
    !> Row value.
    INTEGER, INTENT(OUT) :: index_row
    !> Actual stored value.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: point_value

#include "triplet_includes/GetTriplet.f90"

  END SUBROUTINE GetTripletValues_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compare two triplets based on their index values, first by column and
  !> second by row. Returns A < B.
  PURE FUNCTION CompareTriplets_r(tripA, tripB) RESULT(islessthan)
    !> First triplet.
    TYPE(Triplet_r), INTENT(IN) :: tripA
    !> Second triplet.
    TYPE(Triplet_r), INTENT(IN) :: tripB
    !> A < B.
    LOGICAL :: islessthan

#include "triplet_includes/CompareTriplets.f90"

  END FUNCTION CompareTriplets_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compare two triplets based on their index values (complex), first by
  !> column and second by row. Returns A < B.
  PURE FUNCTION CompareTriplets_c(tripA, tripB) RESULT(islessthan)
    !> First triplet.
    TYPE(Triplet_c), INTENT(IN) :: tripA
    !> Second triplet.
    TYPE(Triplet_c), INTENT(IN) :: tripB
    !> A < B.
    LOGICAL :: islessthan

#include "triplet_includes/CompareTriplets.f90"

  END FUNCTION CompareTriplets_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns an MPI derived data type for a triplet (Real).
  !> We statically store this derived type so that we do not have to recreate
  !> it every time this function is called. Thus this functional call should
  !> add very little overhead.
  FUNCTION GetMPITripletType_r() RESULT(mpi_triplet_type)
    !> MPI Derived Type
    INTEGER :: mpi_triplet_type
    !! Local Data
    INTEGER, DIMENSION(3) :: triplet_sub_types
    INTEGER(KIND = MPI_ADDRESS_KIND), DIMENSION(3) :: triplet_displacement
    INTEGER, DIMENSION(3) :: triplet_block_length
    INTEGER :: bytes_per_int, bytes_per_double
    INTEGER :: ierr

    CALL MPI_Type_extent(MPINTINTEGER, bytes_per_int, ierr)
    CALL MPI_Type_extent(MPINTREAL, bytes_per_double, ierr)
    triplet_block_length = [1, 1, 1]
    triplet_displacement = [0, bytes_per_int, 2 * bytes_per_int]
    triplet_sub_types = [MPINTINTEGER, MPINTINTEGER, MPINTREAL]

    CALL MPI_Type_create_struct(3, triplet_block_length, &
         & triplet_displacement, triplet_sub_types, mpi_triplet_type, ierr)
    CALL MPI_Type_commit(mpi_triplet_type, ierr)

  END FUNCTION GetMPITripletType_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns an MPI derived data type for a triplet (complex).
  !> We statically store this derived type so that we do not have to recreate
  !> it every time this function is called. Thus this functional call should
  !> add very little overhead.
  FUNCTION GetMPITripletType_c() RESULT(mpi_triplet_type)
    !> MPI Derived Type
    INTEGER :: mpi_triplet_type
    !! Local Data
    INTEGER, DIMENSION(3) :: triplet_sub_types
    INTEGER(KIND = MPI_ADDRESS_KIND), DIMENSION(3) :: triplet_displacement
    INTEGER, DIMENSION(3) :: triplet_block_length
    INTEGER :: bytes_per_int, bytes_per_double
    TYPE(Error_t) :: error_check
    LOGICAL :: error_occured
    INTEGER :: ierr

    CALL MPI_Type_extent(MPINTINTEGER, bytes_per_int, ierr)
    CALL MPI_Type_extent(MPINTCOMPLEX, bytes_per_double, ierr)
    triplet_block_length = [1, 1, 1]
    triplet_displacement = [0, bytes_per_int, 2 * bytes_per_int]
    triplet_sub_types = [MPINTINTEGER, MPINTINTEGER, MPINTCOMPLEX]

    CALL MPI_Type_create_struct(3, triplet_block_length, &
         & triplet_displacement, triplet_sub_types, mpi_triplet_type, ierr)
    CALL MPI_Type_commit(mpi_triplet_type, ierr)

    error_occured = CheckMPIError(error_check, "Creation of MPINTCOMPLEX", &
         & ierr, .TRUE.)

  END FUNCTION GetMPITripletType_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a complex triplet to a real triplet.
  SUBROUTINE ConvertTripletToReal(cin_triplet, rout_triplet)
    !> The starting triplet
    TYPE(Triplet_c), INTENT(IN) :: cin_triplet
    !> Real valued triplet.
    TYPE(Triplet_r), INTENT(INOUT) :: rout_triplet

    rout_triplet%index_row = cin_triplet%index_row
    rout_triplet%index_column = cin_triplet%index_column
    rout_triplet%point_value = REAL(cin_triplet%point_value, KIND = NTREAL)
  END SUBROUTINE ConvertTripletToReal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a real triplet to a complex triplet.
  SUBROUTINE ConvertTripletToComplex(rin_triplet, cout_triplet)
    !> The starting triplet.
    TYPE(Triplet_r), INTENT(IN) :: rin_triplet
    !> Complex valued triplet.
    TYPE(Triplet_c), INTENT(INOUT) :: cout_triplet

    cout_triplet%index_row = rin_triplet%index_row
    cout_triplet%index_column = rin_triplet%index_column
    cout_triplet%point_value = CMPLX(rin_triplet%point_value, 0.0_NTREAL, &
         & KIND = NTCOMPLEX)
  END SUBROUTINE ConvertTripletToComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletModule
