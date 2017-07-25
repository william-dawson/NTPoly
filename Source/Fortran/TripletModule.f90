!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Double.
!! Contains a class Triplet and a means to compare them.
MODULE TripletModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL
  USE ErrorModule
  USE ISO_C_BINDING, ONLY : c_int
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a triplet of integer, integer, double.
  !! As this is related to matrix multiplication, the referencing indices are
  !! rows and columns.
  TYPE, PUBLIC :: Triplet_t
     INTEGER(kind=c_int)    :: index_column !< column value.
     INTEGER(kind=c_int)    :: index_row    !< row value.
     REAL(NTREAL) :: point_value  !< actual value at those indices.
  END TYPE Triplet_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetTriplet
  PUBLIC :: GetTripletValues
  PUBLIC :: CompareTriplets
  PUBLIC :: GetMPITripletType
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Flag about whether we've registered an mpi triplet type before
  LOGICAL :: set_mpi_triplet_type = .FALSE.
  !> A derived data type for mpi triplets
  INTEGER, SAVE :: mpi_triplet_type
  !> For error handling in this module
  !type(Error_t) :: error_handler
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  !! @param[inout] this the triplet to set the values of.
  !! @param[in] index_column the column value.
  !! @param[in] index_row the row value.
  !! @param[in] point_value the value at that point.
  PURE SUBROUTINE SetTriplet(this,index_column,index_row,point_value)
    !! Parameters
    TYPE(Triplet_t), INTENT(inout)       :: this
    INTEGER, INTENT(in)                  :: index_column
    INTEGER, INTENT(in)                  :: index_row
    REAL(NTREAL), INTENT(in)   :: point_value

    this%index_column = index_column
    this%index_row    = index_row
    this%point_value  = point_value
  END SUBROUTINE SetTriplet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  !! @param[in] this the triplet to extract the values of.
  !! @param[out] index_column column value.
  !! @param[out] index_row row value.
  !! @param[out] point_value actual stored value.
  PURE SUBROUTINE GetTripletValues(this,index_column,index_row,point_value)
    !! Parameters
    TYPE(Triplet_t), INTENT(in)           :: this
    INTEGER, INTENT(out)                  :: index_column
    INTEGER, INTENT(out)                  :: index_row
    REAL(NTREAL), INTENT(out)   :: point_value

    index_column = this%index_column
    index_row    = this%index_row
    point_value  = this%point_value
  END SUBROUTINE GetTripletValues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compare two triplets based on their index values. Returns A < B.
  !! Used for sorting.
  !! @param[in] tripA first triplet.
  !! @param[in] tripB second triplet.
  !! @return A < B.
  PURE FUNCTION CompareTriplets(tripA, tripB) RESULT(islessthan)
    !! Parameters
    TYPE(Triplet_t), INTENT(in) :: tripA, tripB
    LOGICAL :: islessthan

    IF (tripA%index_column .GT. tripB%index_column) THEN
       islessthan = .TRUE.
    ELSE IF ((tripA%index_column .EQ. tripB%index_column) .AND. &
         & (tripA%index_row .GT. tripB%index_row)) THEN
       islessthan = .TRUE.
    ELSE
       islessthan = .FALSE.
    END IF
  END FUNCTION CompareTriplets
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns an MPI Derived Data Type For A Triplet.
  !! We statically store this derived type so that we don't have to recreate
  !! it every time this function is called. Thus this functional call should
  !! add very little overhead.
  !! @return the derived type
  FUNCTION GetMPITripletType() RESULT(return_mpi_triplet_type)
    !! Parameters
    INTEGER :: return_mpi_triplet_type
    !! Local Data
    INTEGER, DIMENSION(3) :: triplet_sub_types
    INTEGER, DIMENSION(3) :: triplet_displacement
    INTEGER, DIMENSION(3) :: triplet_block_length
    INTEGER :: bytes_per_int
    INTEGER :: bytes_per_double
    INTEGER :: ierr
    !  logical :: error_occurred

    !! Check if already created this data type
    IF (set_mpi_triplet_type) THEN
       return_mpi_triplet_type = mpi_triplet_type
    ELSE
       !! Otherwise make it
       CALL MPI_Type_extent(MPI_INT,bytes_per_int,ierr)
       !error_occurred = CheckAllocError(error_handler, "Creating MPI triplet", &
       !  & ierr, immediate_cleanup_in=.true.)
       CALL MPI_Type_extent(MPINTREAL,bytes_per_double,ierr)
       !error_occurred = CheckAllocError(error_handler, "Creating MPI triplet", &
       !  & ierr, immediate_cleanup_in=.true.)
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
       return_mpi_triplet_type = mpi_triplet_type
    END IF
  END FUNCTION GetMPITripletType
END MODULE TripletModule
