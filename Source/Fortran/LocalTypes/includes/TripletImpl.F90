!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  !! @param[inout] this the triplet to set the values of.
  !! @param[in] index_column the column value.
  !! @param[in] index_row the row value.
  !! @param[in] point_value the value at that point.
  PURE SUBROUTINE SetTriplet(this,index_column,index_row,point_value)
    !! Parameters
    TYPE(INNERTYPE), INTENT(INOUT)       :: this
    INTEGER, INTENT(IN)                  :: index_column
    INTEGER, INTENT(IN)                  :: index_row
    DATATYPE, INTENT(IN)   :: point_value

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
    TYPE(INNERTYPE), INTENT(IN)           :: this
    INTEGER, INTENT(OUT)                  :: index_column
    INTEGER, INTENT(OUT)                  :: index_row
    DATATYPE, INTENT(OUT)   :: point_value

    index_column = this%index_column
    index_row    = this%index_row
    point_value  = this%point_value
  END SUBROUTINE GetTripletValues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compare two triplets based on their index values, first by column and
  !! second by row. Returns A < B.
  !! Used for sorting.
  !! @param[in] tripA first triplet.
  !! @param[in] tripB second triplet.
  !! @return A < B.
  PURE FUNCTION CompareTriplets(tripA, tripB) RESULT(islessthan)
    !! Parameters
    TYPE(INNERTYPE), INTENT(IN) :: tripA, tripB
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
  FUNCTION GetMPITripletType() RESULT(mpi_triplet_type)
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
    CALL MPI_Type_extent(MPIDATATYPE,bytes_per_double,ierr)
    triplet_block_length(1) = 1
    triplet_block_length(2) = 1
    triplet_block_length(3) = 1
    triplet_displacement(1) = 0
    triplet_displacement(2) = bytes_per_int + triplet_displacement(1)
    triplet_displacement(3) = bytes_per_int + triplet_displacement(2)
    triplet_sub_types(1) = MPI_INT
    triplet_sub_types(2) = MPI_INT
    triplet_sub_types(3) = MPIDATATYPE

    CALL MPI_Type_struct(3,triplet_block_length,triplet_displacement,&
         & triplet_sub_types,mpi_triplet_type,ierr)
    CALL MPI_Type_commit(mpi_triplet_type,ierr)
    
  END FUNCTION GetMPITripletType
