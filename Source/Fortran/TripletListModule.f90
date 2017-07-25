!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Lists of triplets.
!! Contains both a methods for sorting lists.
MODULE TripletListModule
  USE MatrixMarketModule
  USE TripletModule, ONLY : Triplet_t, CompareTriplets
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  !! As this is related to matrix multiplication, the referencing indices are
  !! rows and columns.
  TYPE, PUBLIC :: TripletList_t
     !> Internal representation of the data.
     TYPE(Triplet_t), DIMENSION(:), ALLOCATABLE :: DATA
     INTEGER :: CurrentSize
  END TYPE TripletList_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructTripletList
  PUBLIC :: DestructTripletList
  PUBLIC :: ResizeTripletList
  PUBLIC :: AppendToTripletList
  PUBLIC :: SetTripletAt
  PUBLIC :: GetTripletAt
  PUBLIC :: SortTripletList
  PUBLIC :: SymmetrizeTripletList
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  !! @param[inout] this the triplet list to construct.
  !! @param[in] size the length of the triplet list.
  PURE SUBROUTINE ConstructTripletList(this,size_in)
    !! Parameters
    TYPE(TripletList_t), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in
    !! Local data
    INTEGER :: size

    IF (PRESENT(size_in)) THEN
       size = size_in
    ELSE
       size = 0
    END IF

    IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
    this%CurrentSize = size

    ALLOCATE(this%data(size))
  END SUBROUTINE ConstructTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] this the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList(this)
    !! Parameters
    TYPE(TripletList_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
  END SUBROUTINE DestructTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  !! @param[inout] this the triplet list to resize.
  !! @param[in] size.
  PURE SUBROUTINE ResizeTripletList(this, size)
    !! Parameters
    TYPE(TripletList_t), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_t), DIMENSION(:), ALLOCATABLE :: temporary_data

    !! Temporary copy
    ALLOCATE(temporary_data(this%CurrentSize))
    temporary_data = this%DATA(:this%CurrentSize)

    !! Create new memory
    IF (ALLOCATED(this%DATA)) DEALLOCATE(this%DATA)
    ALLOCATE(this%DATA(size))

    !! Copy back
    this%DATA(:this%CurrentSize) = temporary_data

    !! Cleanup
    DEALLOCATE(temporary_data)
  END SUBROUTINE ResizeTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Add a value to the end of the triplet list.
  !! @param[inout] this the triplet list to append to.
  !! @param[in] triplet_value the value to append.
  PURE SUBROUTINE AppendToTripletList(this, triplet_value)
    !! Parameters
    TYPE(TripletList_t), INTENT(inout) :: this
    TYPE(Triplet_t), INTENT(in)        :: triplet_value
    !! Local data
    INTEGER :: new_size

    !! First, check if we need to allocate more memory
    IF (this%CurrentSize+1 .GT. SIZE(this%DATA)) THEN
       IF (SIZE(this%data) .EQ. 0) THEN
          new_size = 1
       ELSE IF (SIZE(this%data) .EQ. 1) THEN
          new_size = 2
       ELSE
          new_size = INT(SIZE(this%data)*1.5)
       END IF
       CALL ResizeTripletList(this,new_size)
    END IF

    !! Append
    this%CurrentSize = this%CurrentSize+1
    this%DATA(this%CurrentSize) = triplet_value

  END SUBROUTINE AppendToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  !! @param[inout] this the triplet list to set.
  !! @param[in] index the index at which to set the triplet.
  !! @param[in] triplet_value the value of the triplet to set.
  PURE SUBROUTINE SetTripletAt(this,index,triplet_value)
    !! Parameters
    TYPE(TripletList_t), INTENT(inout) :: this
    INTEGER(KIND=c_int), INTENT(in)    :: index
    TYPE(Triplet_t), INTENT(in)        :: triplet_value

    this%data(index) = triplet_value
  END SUBROUTINE SetTripletAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  !! @param[in] this the triplet list to get teh value from.
  !! @param[in] index the index from which to get the triplet.
  !! @param[out] triplet_value the extracted triplet value.
  PURE SUBROUTINE GetTripletAt(this,index,triplet_value)
    !! Parameters
    TYPE(TripletList_t), INTENT(in) :: this
    INTEGER(kind=c_int), INTENT(in) :: index
    TYPE(Triplet_t), INTENT(out)    :: triplet_value

    triplet_value = this%data(index)
  END SUBROUTINE GetTripletAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !! Implementation is based on bucket sort. This is why it needs the number of
  !! matrix columns. Bubble sort is used within a bucket.
  !! @param[in] input_list list to be sorted.
  !! @param[in] matrix_columns this is the highest column value in the list.
  !! @param[in] bubble_in false if you don't need the final bubble sort.
  !! @param[out] sorted_list a now sorted version of the list. This routine
  !! will allocate it.
  PURE SUBROUTINE SortTripletList(input_list, matrix_columns, sorted_list, &
       & bubble_in)
    !! Parameters
    TYPE(TripletList_t), INTENT(in)  :: input_list
    INTEGER, INTENT(in) :: matrix_columns
    TYPE(TripletList_t), INTENT(out) :: sorted_list
    LOGICAL, OPTIONAL, INTENT(in) :: bubble_in
    !! Local Data
    LOGICAL :: bubble
    TYPE(Triplet_t) :: temporary
    LOGICAL :: swap_occured
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
    !! Counters and temporary variables
    INTEGER :: counter
    INTEGER :: temp_index
    INTEGER :: alloc_stat
    INTEGER :: list_length

    IF (PRESENT(bubble_in)) THEN
       bubble = bubble_in
    ELSE
       bubble = .TRUE.
    END IF

    !! Data Allocation
    !list_length = SIZE(input_list%data)
    list_length = input_list%CurrentSize
    CALL ConstructTripletList(sorted_list,list_length)
    !allocate(sorted_list%(size(input_list%data)), stat=alloc_stat)
    ALLOCATE(values_per_row(matrix_columns), stat=alloc_stat)
    ALLOCATE(offset_array(matrix_columns), stat=alloc_stat)
    ALLOCATE(inserted_per_row(matrix_columns), stat=alloc_stat)

    !! Initial one dimensional sort
    values_per_row = 0
    inserted_per_row = 0

    !! Do a first pass bucket sort
    !DO counter = LBOUND(input_list%data(:),dim=1), UBOUND(input_list%data,dim=1)
    DO counter = 1, input_list%CurrentSize
       values_per_row(input_list%data(counter)%index_column) = &
            & values_per_row(input_list%data(counter)%index_column) + 1
    END DO
    offset_array(1) = 1
    DO counter = 2, UBOUND(offset_array,dim=1)
       offset_array(counter) = offset_array(counter-1) + &
            & values_per_row(counter-1)
    END DO
    !DO counter = LBOUND(input_list%data,dim=1), UBOUND(input_list%data,dim=1)
    DO counter = 1, input_list%CurrentSize
       temp_index = input_list%data(counter)%index_column
       sorted_list%data(offset_array(temp_index)+inserted_per_row(temp_index))=&
            & input_list%data(counter)
       inserted_per_row(temp_index) = inserted_per_row(temp_index) + 1
    END DO

    !! Finish with bubble sort
    !! Not necessary for transpoing.
    swap_occured = .TRUE.
    IF (bubble) THEN
       DO WHILE (swap_occured .EQV. .TRUE.)
          swap_occured = .FALSE.
          DO counter = 2, sorted_list%CurrentSize
             IF (CompareTriplets(sorted_list%data(counter-1), &
                  & sorted_list%data(counter))) THEN
                temporary = sorted_list%data(counter)
                sorted_list%data(counter) = sorted_list%data(counter-1)
                sorted_list%data(counter-1) = temporary
                swap_occured = .TRUE.
             END IF
          END DO
       END DO
    END IF

    !! Cleanup
    DEALLOCATE(values_per_row)
    DEALLOCATE(offset_array)
    DEALLOCATE(inserted_per_row)
  END SUBROUTINE SortTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] input_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList(triplet_list, pattern_type)
    !! Parameters
    TYPE(TripletList_t), INTENT(INOUT)  :: triplet_list
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(Triplet_t) :: temporary, temporary_transpose
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
END MODULE TripletListModule
