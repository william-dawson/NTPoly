!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Value.
MODULE TripletModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE ISO_C_BINDING, ONLY : c_int
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base type for triplets, which has no data value.
  TYPE, ABSTRACT, PUBLIC :: Triplet
    INTEGER(kind=c_int) :: index_column !< column value.
    INTEGER(kind=c_int) :: index_row    !< row value.
  CONTAINS
    PROCEDURE, PASS :: LessThan => Compare_base
    PROCEDURE :: GetMPIType => GetMPIType_base
    PROCEDURE :: Transpose => Transpose_base
    PROCEDURE :: Conjg => Conjg_base
    PROCEDURE :: Scale => Scale_base
    PROCEDURE :: ReadFromFile => Read_base
    PROCEDURE :: WriteToFile => Write_base
  END TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a triplet of integer, integer, double.
  TYPE, EXTENDS(Triplet), PUBLIC :: Triplet_r
     REAL(NTREAL) :: point_value !< actual value at those indices.
  CONTAINS
     PROCEDURE :: SetTriplet => SetTriplet_r
     PROCEDURE :: GetTriplet => GetTriplet_r
     PROCEDURE :: ReadFromFile => Read_r
     PROCEDURE :: WriteToFile => Write_r
  END TYPE Triplet_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a triplet of integer, integer, complex.
  TYPE, EXTENDS(Triplet), PUBLIC :: Triplet_c
     COMPLEX(NTCOMPLEX) :: point_value !< actual value at those indices.
  CONTAINS
     PROCEDURE :: SetTriplet => SetTriplet_c
     PROCEDURE :: GetTriplet => GetTriplet_c
     PROCEDURE :: ReadFromFile => Read_c
     PROCEDURE :: WriteToFile => Write_c
  END TYPE Triplet_c
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns an MPI Derived Data Type For A Triplet.
  !! We statically store this derived type so that we don't have to recreate
  !! it every time this function is called. Thus this functional call should
  !! add very little overhead.
  !! @return the derived type
  FUNCTION GetMPIType_base(this) RESULT(mpi_triplet_type)
    !! Parameters
    CLASS(Triplet), INTENT(IN) :: this
    INTEGER :: mpi_triplet_type
    !! Local Data
    INTEGER, DIMENSION(3) :: triplet_sub_types
    INTEGER, DIMENSION(3) :: triplet_displacement
    INTEGER, DIMENSION(3) :: triplet_block_length
    INTEGER :: bytes_per_int
    INTEGER :: ierr

    CALL MPI_Type_extent(MPI_INT,bytes_per_int,ierr)
    triplet_block_length(1) = 1
    triplet_block_length(2) = 1
    triplet_block_length(3) = 1
    triplet_displacement(1) = 0
    triplet_displacement(2) = bytes_per_int + triplet_displacement(1)
    triplet_displacement(3) = bytes_per_int + triplet_displacement(2)
    triplet_sub_types(1) = MPI_INT
    triplet_sub_types(2) = MPI_INT

    SELECT TYPE(this)
    CLASS IS(Triplet_r)
       triplet_sub_types(3) = MPINTREAL
    CLASS IS(Triplet_c)
       triplet_sub_types(3) = MPINTCOMPLEX
    END SELECT

    CALL MPI_Type_struct(3,triplet_block_length,triplet_displacement,&
         & triplet_sub_types,mpi_triplet_type,ierr)
    CALL MPI_Type_commit(mpi_triplet_type,ierr)

  END FUNCTION GetMPIType_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compare two triplets based on their index values, first by column and
  !! second by row. Returns A < B.
  !! Used for sorting.
  !! @param[in] tripA first triplet.
  !! @param[in] tripB second triplet.
  !! @return A < B.
  PURE FUNCTION Compare_base(tripA, tripB) RESULT(islessthan)
    !! Parameters
    CLASS(Triplet), INTENT(IN) :: tripA, tripB
    LOGICAL :: islessthan

    IF (tripA%index_column .GT. tripB%index_column) THEN
       islessthan = .TRUE.
    ELSE IF ((tripA%index_column .EQ. tripB%index_column) .AND. &
         & (tripA%index_row .GT. tripB%index_row)) THEN
       islessthan = .TRUE.
    ELSE
       islessthan = .FALSE.
    END IF
  END FUNCTION Compare_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Exchange the rows and columns of a triplet.
  !! @param[inout] this the triplet to transpose.
  PURE SUBROUTINE Transpose_base(this)
    !! Parameters
    CLASS(Triplet), INTENT(INOUT) :: this
    INTEGER :: temp

    temp = this%index_row
    this%index_row =this%index_column
    this%index_column = temp

  END SUBROUTINE Transpose_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read a line of an open file handle into this triplet.
  !! @param[inout] this the triplet to read in
  !! @param[inout] file_handle an open file handle sitting at a line to read.
  SUBROUTINE Read_base(this, file_handle)
    !! Parameters
    CLASS(Triplet), INTENT(INOUT) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    READ(file_handle,*) this%index_row, this%index_column

  END SUBROUTINE Read_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read a line of an open file handle into this triplet.
  !! @param[inout] this the triplet to read in
  !! @param[inout] file_handle an open file handle sitting at a line to read.
  SUBROUTINE Read_r(this, file_handle)
    !! Parameters
    CLASS(Triplet_r), INTENT(INOUT) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    READ(file_handle,*) this%index_row, this%index_column, &
         & this%point_value

  END SUBROUTINE Read_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read a line of an open file handle into this triplet.
  !! @param[inout] this the triplet to read in
  !! @param[inout] file_handle an open file handle sitting at a line to read.
  SUBROUTINE Read_c(this, file_handle)
    !! Parameters
    CLASS(Triplet_c), INTENT(INOUT) :: this
    INTEGER, INTENT(INOUT) :: file_handle
    !! Local
    REAL(NTREAL) :: real_val, comp_val

    READ(file_handle,*) this%index_row, this%index_column, real_val, comp_val
    this%point_value = CMPLX(real_val, comp_val, KIND=NTCOMPLEX)

  END SUBROUTINE Read_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read a line of an open file handle into this triplet.
  !! @param[inout] this the triplet to read in
  !! @param[inout] file_handle an open file handle sitting at a line to read.
  SUBROUTINE Write_base(this, file_handle)
    !! Parameters
    CLASS(Triplet), INTENT(IN) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    WRITE(file_handle,*) this%index_row, this%index_column

  END SUBROUTINE Write_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read a line of an open file handle into this triplet.
  !! @param[inout] this the triplet to read in
  !! @param[inout] file_handle an open file handle sitting at a line to read.
  SUBROUTINE Write_r(this, file_handle)
    !! Parameters
    CLASS(Triplet_r), INTENT(IN) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    WRITE(file_handle,*) this%index_row, this%index_column, this%point_value

  END SUBROUTINE Write_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read a line of an open file handle into this triplet.
  !! @param[inout] this the triplet to read in
  !! @param[inout] file_handle an open file handle sitting at a line to read.
  SUBROUTINE Write_c(this, file_handle)
    !! Parameters
    CLASS(Triplet_c), INTENT(IN) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    WRITE(file_handle,*) this%index_row, this%index_column, &
         & REAL(this%point_value), AIMAG(this%point_value)

  END SUBROUTINE Write_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  !! @param[inout] this the triplet to set the values of.
  !! @param[in] index_column the column value.
  !! @param[in] index_row the row value.
  !! @param[in] point_value the value at that point.
  PURE SUBROUTINE SetTriplet_r(this,index_column,index_row,point_value)
    !! Parameters
    CLASS(Triplet_r), INTENT(INOUT) :: this
    INTEGER, INTENT(IN)             :: index_column
    INTEGER, INTENT(IN)             :: index_row
    REAL(NTREAL), INTENT(IN)        :: point_value

    this%index_column = index_column
    this%index_row    = index_row
    this%point_value  = point_value

  END SUBROUTINE SetTriplet_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  !! @param[in] this the triplet to extract the values of.
  !! @param[out] index_column column value.
  !! @param[out] index_row row value.
  !! @param[out] point_value actual stored value.
  PURE SUBROUTINE GetTriplet_r(this,index_column,index_row,point_value)
    !! Parameters
    CLASS(Triplet_r), INTENT(IN) :: this
    INTEGER, INTENT(OUT)        :: index_column
    INTEGER, INTENT(OUT)        :: index_row
    REAL(NTREAL), INTENT(OUT)   :: point_value

    index_column = this%index_column
    index_row    = this%index_row
    point_value  = this%point_value

  END SUBROUTINE GetTriplet_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  !! @param[inout] this the triplet to set the values of.
  !! @param[in] index_column the column value.
  !! @param[in] index_row the row value.
  !! @param[in] point_value the value at that point.
  PURE SUBROUTINE SetTriplet_c(this,index_column,index_row,point_value)
    !! Parameters
    CLASS(Triplet_c), INTENT(INOUT) :: this
    INTEGER, INTENT(IN)             :: index_column
    INTEGER, INTENT(IN)             :: index_row
    COMPLEX(NTCOMPLEX), INTENT(IN)  :: point_value

    this%index_column = index_column
    this%index_row    = index_row
    this%point_value  = point_value

  END SUBROUTINE SetTriplet_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  !! @param[in] this the triplet to extract the values of.
  !! @param[out] index_column column value.
  !! @param[out] index_row row value.
  !! @param[out] point_value actual stored value.
  PURE SUBROUTINE GetTriplet_c(this,index_column,index_row,point_value)
    !! Parameters
    CLASS(Triplet_c), INTENT(IN)    :: this
    INTEGER, INTENT(OUT)            :: index_column
    INTEGER, INTENT(OUT)            :: index_row
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: point_value

    index_column = this%index_column
    index_row    = this%index_row
    point_value  = this%point_value

  END SUBROUTINE GetTriplet_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE Conjg_base(this)
    !! Parameters
    CLASS(Triplet), INTENT(INOUT) :: this

    SELECT TYPE(this)
    CLASS IS (Triplet_c)
      this%point_value = CONJG(this%point_value)
    END SELECT
  END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE Scale_base(this, alpha)
    !! Parameters
    CLASS(Triplet), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: alpha

    SELECT TYPE(this)
    CLASS IS (Triplet_r)
      this%point_value = alpha * this%point_value
    CLASS IS (Triplet_c)
      this%point_value = alpha * this%point_value
    END SELECT
  END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletModule
