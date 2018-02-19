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
     TYPE(Triplet_t), POINTER :: DATA
  END TYPE Triplet_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetTriplet_wrp
  PUBLIC :: GetTripletValues_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  PURE SUBROUTINE SetTriplet_wrp(ih_this,index_column,index_row,point_value) &
       & bind(c,name="SetTriplet_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    REAL(NTREAL), INTENT(IN) :: point_value
    TYPE(Triplet_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetTriplet(h_this%data,index_column,index_row,point_value)
  END SUBROUTINE SetTriplet_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  PURE SUBROUTINE GetTripletValues_wrp(ih_this,index_column,index_row, &
       & point_value) bind(c,name="GetTripletValues_wrp")
    INTEGER(kind=c_int), INTENT(IN)     :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)    :: index_column
    INTEGER(kind=c_int), INTENT(OUT)    :: index_row
    REAL(NTREAL), INTENT(OUT) :: point_value
    TYPE(Triplet_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL GetTripletValues(h_this%data,index_column,index_row,point_value)
  END SUBROUTINE GetTripletValues_wrp
END MODULE TripletModule_wrp
