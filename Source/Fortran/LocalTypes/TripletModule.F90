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
  !> Set the values of a triplet.
  !! @param[inout] this the triplet to set the values of.
  !! @param[in] index_column the column value.
  !! @param[in] index_row the row value.
  !! @param[in] point_value the value at that point.
  PURE SUBROUTINE SetTriplet_r(this,index_column,index_row,point_value)
    !! Parameters
    TYPE(Triplet_r), INTENT(INOUT)       :: this
    INTEGER, INTENT(IN)                  :: index_column
    INTEGER, INTENT(IN)                  :: index_row
    REAL(NTREAL), INTENT(IN)   :: point_value

    INCLUDE "triplet_includes/SetTriplet.f90"

  END SUBROUTINE SetTriplet_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  !! @param[in] this the triplet to extract the values of.
  !! @param[out] index_column column value.
  !! @param[out] index_row row value.
  !! @param[out] point_value actual stored value.
  PURE SUBROUTINE GetTripletValues_r(this,index_column,index_row,point_value)
    !! Parameters
    TYPE(Triplet_r), INTENT(IN)           :: this
    INTEGER, INTENT(OUT)                  :: index_column
    INTEGER, INTENT(OUT)                  :: index_row
    REAL(NTREAL), INTENT(OUT)   :: point_value

    INCLUDE "triplet_includes/GetTriplet.f90"

  END SUBROUTINE GetTripletValues_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compare two triplets based on their index values, first by column and
  !! second by row. Returns A < B.
  !! Used for sorting.
  !! @param[in] tripA first triplet.
  !! @param[in] tripB second triplet.
  !! @return A < B.
  PURE FUNCTION CompareTriplets_r(tripA, tripB) RESULT(islessthan)
    !! Parameters
    TYPE(Triplet_r), INTENT(IN) :: tripA, tripB
    LOGICAL :: islessthan

    INCLUDE "triplet_includes/CompareTriplets.f90"

  END FUNCTION CompareTriplets_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns an MPI Derived Data Type For A Triplet.
  !! We statically store this derived type so that we don't have to recreate
  !! it every time this function is called. Thus this functional call should
  !! add very little overhead.
  !! @return the derived type
  FUNCTION GetMPITripletType_r() RESULT(mpi_triplet_type)
    !! Parameters
    INTEGER :: mpi_triplet_type
    !! Local Data
    INTEGER, DIMENSION(3) :: triplet_sub_types
    INTEGER, DIMENSION(3) :: triplet_displacement
    INTEGER, DIMENSION(3) :: triplet_block_length
    INTEGER :: bytes_per_int
    INTEGER :: bytes_per_double
    INTEGER :: ierr

    CALL MPI_Type_extent(MPI_INT,bytes_per_int,ierr)
    CALL MPI_Type_extent(MPINTREAL,bytes_per_double,ierr)
    triplet_block_length(1) = 1
    triplet_block_length(2) = 1
    triplet_block_length(3) = 1
    triplet_displacement(1) = 0
    triplet_displacement(2) = bytes_per_int + triplet_displacement(1)
    triplet_displacement(3) = bytes_per_int + triplet_displacement(2)
    triplet_sub_types(1) = MPI_INT
    triplet_sub_types(2) = MPI_INT
    triplet_sub_types(3) = MPINTREAL

    CALL MPI_Type_struct(3,triplet_block_length,triplet_displacement,&
         & triplet_sub_types,mpi_triplet_type,ierr)
    CALL MPI_Type_commit(mpi_triplet_type,ierr)

  END FUNCTION GetMPITripletType_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  !! @param[inout] this the triplet to set the values of.
  !! @param[in] index_column the column value.
  !! @param[in] index_row the row value.
  !! @param[in] point_value the value at that point.
  PURE SUBROUTINE SetTriplet_c(this,index_column,index_row,point_value)
    !! Parameters
    TYPE(Triplet_c), INTENT(INOUT)       :: this
    INTEGER, INTENT(IN)                  :: index_column
    INTEGER, INTENT(IN)                  :: index_row
    COMPLEX(NTCOMPLEX), INTENT(IN)   :: point_value

    INCLUDE "triplet_includes/SetTriplet.f90"

  END SUBROUTINE SetTriplet_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  !! @param[in] this the triplet to extract the values of.
  !! @param[out] index_column column value.
  !! @param[out] index_row row value.
  !! @param[out] point_value actual stored value.
  PURE SUBROUTINE GetTripletValues_c(this,index_column,index_row,point_value)
    !! Parameters
    TYPE(Triplet_c), INTENT(IN)           :: this
    INTEGER, INTENT(OUT)                  :: index_column
    INTEGER, INTENT(OUT)                  :: index_row
    COMPLEX(NTCOMPLEX), INTENT(OUT)   :: point_value

    INCLUDE "triplet_includes/GetTriplet.f90"

  END SUBROUTINE GetTripletValues_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compare two triplets based on their index values, first by column and
  !! second by row. Returns A < B.
  !! Used for sorting.
  !! @param[in] tripA first triplet.
  !! @param[in] tripB second triplet.
  !! @return A < B.
  PURE FUNCTION CompareTriplets_c(tripA, tripB) RESULT(islessthan)
    !! Parameters
    TYPE(Triplet_c), INTENT(IN) :: tripA, tripB
    LOGICAL :: islessthan

    INCLUDE "triplet_includes/CompareTriplets.f90"

  END FUNCTION CompareTriplets_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns an MPI Derived Data Type For A Triplet.
  !! We statically store this derived type so that we don't have to recreate
  !! it every time this function is called. Thus this functional call should
  !! add very little overhead.
  !! @return the derived type
  FUNCTION GetMPITripletType_c() RESULT(mpi_triplet_type)
    !! Parameters
    INTEGER :: mpi_triplet_type
    !! Local Data
    INTEGER, DIMENSION(3) :: triplet_sub_types
    INTEGER, DIMENSION(3) :: triplet_displacement
    INTEGER, DIMENSION(3) :: triplet_block_length
    INTEGER :: bytes_per_int
    INTEGER :: bytes_per_double
    INTEGER :: ierr

    CALL MPI_Type_extent(MPI_INT,bytes_per_int,ierr)
    CALL MPI_Type_extent(MPINTCOMPLEX,bytes_per_double,ierr)
    triplet_block_length(1) = 1
    triplet_block_length(2) = 1
    triplet_block_length(3) = 1
    triplet_displacement(1) = 0
    triplet_displacement(2) = bytes_per_int + triplet_displacement(1)
    triplet_displacement(3) = bytes_per_int + triplet_displacement(2)
    triplet_sub_types(1) = MPI_INT
    triplet_sub_types(2) = MPI_INT
    triplet_sub_types(3) = MPINTCOMPLEX

    CALL MPI_Type_struct(3,triplet_block_length,triplet_displacement,&
         & triplet_sub_types,mpi_triplet_type,ierr)
    CALL MPI_Type_commit(mpi_triplet_type,ierr)

  END FUNCTION GetMPITripletType_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletModule
