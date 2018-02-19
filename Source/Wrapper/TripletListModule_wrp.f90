!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the triplet list data type.
MODULE TripletListModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE TripletListModule, ONLY : TripletList_t, ConstructTripletList, &
       & AppendToTripletList, ResizeTripletList, &
       & DestructTripletList, SetTripletAt, GetTripletAt, SortTripletList, &
       & GetTripletListSize
  USE TripletModule, ONLY : Triplet_t
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the triplet list data type.
  TYPE, PUBLIC :: TripletList_wrp
     !> Actual data.
     TYPE(TripletList_t), POINTER :: DATA
  END TYPE TripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructTripletList_wrp
  PUBLIC :: DestructTripletList_wrp
  PUBLIC :: ResizeTripletList_wrp
  PUBLIC :: AppendToTripletList_wrp
  PUBLIC :: SetTripletAt_wrp
  PUBLIC :: GetTripletAt_wrp
  PUBLIC :: SortTripletList_wrp
  PUBLIC :: GetTripletListSize_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the triplet list constructor.
  PURE SUBROUTINE ConstructTripletList_wrp(ih_this, size) &
       & bind(c,name="ConstructTripletList_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: size
    TYPE(TripletList_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructTripletList(h_this%data,size)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] ih_this handle to the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList_wrp(ih_this) &
       & bind(c,name="DestructTripletList_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(TripletList_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructTripletList(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  PURE SUBROUTINE ResizeTripletList_wrp(ih_this, size) &
       & bind(c,name="ResizeTripletList_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: size
    TYPE(TripletList_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ResizeTripletList(h_this%data,size)
  END SUBROUTINE ResizeTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  PURE SUBROUTINE AppendToTripletList_wrp(ih_this, index_column, index_row, &
       & point_value) bind(c,name="AppendToTripletList_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    REAL(NTREAL), INTENT(IN) :: point_value
    TYPE(TripletList_wrp) :: h_this
    TYPE(Triplet_t) :: temp

    temp%index_column = index_column
    temp%index_row = index_row
    temp%point_value = point_value

    h_this = TRANSFER(ih_this,h_this)
    CALL AppendToTripletList(h_this%data,temp)
  END SUBROUTINE AppendToTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  PURE SUBROUTINE SetTripletAt_wrp(ih_this, index, index_column, index_row, &
       & point_value) bind(c,name="SetTripletAt_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    REAL(NTREAL), INTENT(IN) :: point_value
    TYPE(TripletList_wrp) :: h_this
    TYPE(Triplet_t) :: temp

    temp%index_column = index_column
    temp%index_row = index_row
    temp%point_value = point_value

    h_this = TRANSFER(ih_this,h_this)
    CALL SetTripletAt(h_this%data,index,temp)
  END SUBROUTINE SetTripletAt_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  PURE SUBROUTINE GetTripletAt_wrp(ih_this, index, index_column, index_row, &
       & point_value) bind(c,name="GetTripletAt_wrp")
    INTEGER(kind=c_int), INTENT(IN)     :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)    :: index
    INTEGER(kind=c_int), INTENT(OUT)    :: index_column
    INTEGER(kind=c_int), INTENT(OUT)    :: index_row
    REAL(NTREAL), INTENT(OUT) :: point_value
    TYPE(TripletList_wrp) :: h_this
    TYPE(Triplet_t) :: temp

    h_this = TRANSFER(ih_this,h_this)
    CALL GetTripletAt(h_this%data,index,temp)

    index_column = temp%index_column
    index_row = temp%index_row
    point_value = temp%point_value
  END SUBROUTINE GetTripletAt_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  PURE SUBROUTINE SortTripletList_wrp(ih_this, matrix_columns, ih_sorted) &
    & bind(c,name="SortTripletList_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: matrix_columns
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_sorted(SIZE_wrp)
    TYPE(TripletList_wrp) :: h_this
    TYPE(TripletList_wrp) :: h_sorted

    h_this = TRANSFER(ih_this,h_this)
    ALLOCATE(h_sorted%data)

    CALL SortTripletList(h_this%data, matrix_columns, h_sorted%data)

    ih_sorted = TRANSFER(h_sorted,ih_sorted)
  END SUBROUTINE SortTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  PURE FUNCTION GetTripletListSize_wrp(ih_this) RESULT(list_size) &
       & bind(c,name="GetTripletListSize_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int) :: list_size
    TYPE(TripletList_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    list_size = GetTripletListSize(h_this%data)
  END FUNCTION GetTripletListSize_wrp
END MODULE TripletListModule_wrp
