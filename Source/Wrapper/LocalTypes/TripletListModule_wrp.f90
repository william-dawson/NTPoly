!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the triplet list data type.
MODULE TripletListModule_wrp
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, SortTripletList
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the triplet list data type.
  TYPE, PUBLIC :: TripletList_r_wrp
     !> Actual data.
     TYPE(TripletList_r), POINTER :: DATA
  END TYPE TripletList_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the triplet list data type.
  TYPE, PUBLIC :: TripletList_c_wrp
     !> Actual data.
     TYPE(TripletList_c), POINTER :: DATA
  END TYPE TripletList_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructTripletList_r_wrp
  PUBLIC :: DestructTripletList_r_wrp
  PUBLIC :: ResizeTripletList_r_wrp
  PUBLIC :: AppendToTripletList_r_wrp
  PUBLIC :: SetTripletAt_r_wrp
  PUBLIC :: GetTripletAt_r_wrp
  PUBLIC :: SortTripletList_r_wrp
  PUBLIC :: GetTripletListSize_r_wrp
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructTripletList_c_wrp
  PUBLIC :: DestructTripletList_c_wrp
  PUBLIC :: ResizeTripletList_c_wrp
  PUBLIC :: AppendToTripletList_c_wrp
  PUBLIC :: SetTripletAt_c_wrp
  PUBLIC :: GetTripletAt_c_wrp
  PUBLIC :: SortTripletList_c_wrp
  PUBLIC :: GetTripletListSize_c_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the triplet list constructor.
  PURE SUBROUTINE ConstructTripletList_r_wrp(ih_this, size) &
       & bind(c,name="ConstructTripletList_r_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: size
    TYPE(TripletList_r_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL h_this%data%Init(size)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructTripletList_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] ih_this handle to the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList_r_wrp(ih_this) &
       & bind(c,name="DestructTripletList_r_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(TripletList_r_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Destruct
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructTripletList_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  PURE SUBROUTINE ResizeTripletList_r_wrp(ih_this, size) &
       & bind(c,name="ResizeTripletList_r_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: size
    TYPE(TripletList_r_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Resize(size)
  END SUBROUTINE ResizeTripletList_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  PURE SUBROUTINE AppendToTripletList_r_wrp(ih_this, index_column, index_row, &
       & point_value) bind(c,name="AppendToTripletList_r_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    REAL(NTREAL), INTENT(IN) :: point_value
    TYPE(TripletList_r_wrp) :: h_this
    TYPE(Triplet_r) :: temp

    temp%index_column = index_column
    temp%index_row = index_row
    temp%point_value = point_value

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Append(temp)
  END SUBROUTINE AppendToTripletList_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  PURE SUBROUTINE SetTripletAt_r_wrp(ih_this, index, index_column, index_row, &
       & point_value) bind(c,name="SetTripletAt_r_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    REAL(NTREAL), INTENT(IN) :: point_value
    TYPE(TripletList_r_wrp) :: h_this
    TYPE(Triplet_r) :: temp

    temp%index_column = index_column
    temp%index_row = index_row
    temp%point_value = point_value

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Set(index, temp)
  END SUBROUTINE SetTripletAt_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  PURE SUBROUTINE GetTripletAt_r_wrp(ih_this, index, index_column, index_row, &
       & point_value) bind(c,name="GetTripletAt_r_wrp")
    INTEGER(kind=c_int), INTENT(IN)     :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)    :: index
    INTEGER(kind=c_int), INTENT(OUT)    :: index_column
    INTEGER(kind=c_int), INTENT(OUT)    :: index_row
    REAL(NTREAL), INTENT(OUT) :: point_value
    TYPE(TripletList_r_wrp) :: h_this
    TYPE(Triplet_r) :: temp

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Get(index, temp)

    index_column = temp%index_column
    index_row = temp%index_row
    point_value = temp%point_value
  END SUBROUTINE GetTripletAt_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  PURE SUBROUTINE SortTripletList_r_wrp(ih_this, matrix_columns, matrix_rows, &
       & ih_sorted) bind(c,name="SortTripletList_r_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: matrix_columns
    INTEGER(kind=c_int), INTENT(IN)    :: matrix_rows
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_sorted(SIZE_wrp)
    TYPE(TripletList_r_wrp) :: h_this
    TYPE(TripletList_r_wrp) :: h_sorted

    h_this = TRANSFER(ih_this,h_this)
    ALLOCATE(h_sorted%data)

    CALL SortTripletList(h_this%data, matrix_columns, matrix_rows, &
         & h_sorted%data)

    ih_sorted = TRANSFER(h_sorted,ih_sorted)
  END SUBROUTINE SortTripletList_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  PURE FUNCTION GetTripletListSize_r_wrp(ih_this) RESULT(list_size) &
       & bind(c,name="GetTripletListSize_r_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int) :: list_size
    TYPE(TripletList_r_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    list_size = h_this%data%GetSize()
  END FUNCTION GetTripletListSize_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the triplet list constructor.
  PURE SUBROUTINE ConstructTripletList_c_wrp(ih_this, size) &
       & bind(c,name="ConstructTripletList_c_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: size
    TYPE(TripletList_c_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL h_this%data%Init(size)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructTripletList_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] ih_this handle to the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList_c_wrp(ih_this) &
       & bind(c,name="DestructTripletList_c_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(TripletList_c_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Destruct
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructTripletList_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  PURE SUBROUTINE ResizeTripletList_c_wrp(ih_this, size) &
       & bind(c,name="ResizeTripletList_c_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: size
    TYPE(TripletList_c_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Resize(size)
  END SUBROUTINE ResizeTripletList_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  PURE SUBROUTINE AppendToTripletList_c_wrp(ih_this, index_column, index_row, &
       & point_value) bind(c,name="AppendToTripletList_c_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    COMPLEX(NTCOMPLEX), INTENT(IN) :: point_value
    TYPE(TripletList_c_wrp) :: h_this
    TYPE(Triplet_c) :: temp

    temp%index_column = index_column
    temp%index_row = index_row
    temp%point_value = point_value

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Append(temp)
  END SUBROUTINE AppendToTripletList_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  PURE SUBROUTINE SetTripletAt_c_wrp(ih_this, index, index_column, index_row, &
       & point_value) bind(c,name="SetTripletAt_c_wrp")
    INTEGER(kind=c_int), INTENT(INOUT)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    COMPLEX(NTCOMPLEX), INTENT(IN) :: point_value
    TYPE(TripletList_c_wrp) :: h_this
    TYPE(Triplet_c) :: temp

    temp%index_column = index_column
    temp%index_row = index_row
    temp%point_value = point_value

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Set(index, temp)
  END SUBROUTINE SetTripletAt_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  PURE SUBROUTINE GetTripletAt_c_wrp(ih_this, index, index_column, index_row, &
       & point_value) bind(c,name="GetTripletAt_c_wrp")
    INTEGER(kind=c_int), INTENT(IN)     :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)    :: index
    INTEGER(kind=c_int), INTENT(OUT)    :: index_column
    INTEGER(kind=c_int), INTENT(OUT)    :: index_row
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: point_value
    TYPE(TripletList_c_wrp) :: h_this
    TYPE(Triplet_c) :: temp

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%Get(index, temp)

    index_column = temp%index_column
    index_row = temp%index_row
    point_value = temp%point_value
  END SUBROUTINE GetTripletAt_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  PURE SUBROUTINE SortTripletList_c_wrp(ih_this, matrix_columns, matrix_rows, &
       & ih_sorted) bind(c,name="SortTripletList_c_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: matrix_columns
    INTEGER(kind=c_int), INTENT(IN)    :: matrix_rows
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_sorted(SIZE_wrp)
    TYPE(TripletList_c_wrp) :: h_this
    TYPE(TripletList_c_wrp) :: h_sorted

    h_this = TRANSFER(ih_this,h_this)
    ALLOCATE(h_sorted%data)

    CALL SortTripletList(h_this%data, matrix_columns, matrix_rows, &
         & h_sorted%data)

    ih_sorted = TRANSFER(h_sorted,ih_sorted)
  END SUBROUTINE SortTripletList_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  PURE FUNCTION GetTripletListSize_c_wrp(ih_this) RESULT(list_size) &
       & bind(c,name="GetTripletListSize_c_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int) :: list_size
    TYPE(TripletList_c_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    list_size = h_this%data%GetSize()
  END FUNCTION GetTripletListSize_c_wrp
END MODULE TripletListModule_wrp
