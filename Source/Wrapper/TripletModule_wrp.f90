!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the triplet module for calling from other languages.
MODULE TripletModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE TripletModule, ONLY: Triplet_t, SetTriplet, GetTripletValues
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the triplet data type.
  TYPE, PUBLIC :: Triplet_wrp
     !> Actual data.
     TYPE(Triplet_t), POINTER :: DATA
  END TYPE Triplet_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetTriplet_wrp
  PUBLIC :: GetTripletValues_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  !! @param[inout] ih_this handle to the triplet to set the values of.
  !! @param[in] index_column the column value.
  !! @param[in] index_row the row value.
  !! @param[in] point_value the value at that point.
  PURE SUBROUTINE SetTriplet_wrp(ih_this,index_column,index_row,point_value) &
       & bind(c,name="SetTriplet_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    REAL(NTREAL), INTENT(IN) :: point_value
    TYPE(Triplet_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetTriplet(h_this%data,index_column,index_row,point_value)
  END SUBROUTINE SetTriplet_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  !! @param[in] ih_this the triplet to extract the values of.
  !! @param[out] index_column column value.
  !! @param[out] index_row row value.
  !! @param[out] point_value actual stored value.
  PURE SUBROUTINE GetTripletValues_wrp(ih_this,index_column,index_row, &
       & point_value) bind(c,name="GetTripletValues_wrp")
    INTEGER(kind=c_int), INTENT(IN)     :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(out)    :: index_column
    INTEGER(kind=c_int), INTENT(out)    :: index_row
    REAL(NTREAL), INTENT(out) :: point_value
    TYPE(Triplet_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL GetTripletValues(h_this%data,index_column,index_row,point_value)
  END SUBROUTINE GetTripletValues_wrp
END MODULE TripletModule_wrp
