!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the triplet module for calling from other languages.
MODULE TripletModule_wrp
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE TripletModule, ONLY: Triplet_r, Triplet_c
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the triplet data type.
  TYPE, PUBLIC :: Triplet_r_wrp
     TYPE(Triplet_r), POINTER :: DATA
  END TYPE Triplet_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the triplet data type.
  TYPE, PUBLIC :: Triplet_c_wrp
     TYPE(Triplet_c), POINTER :: DATA
  END TYPE Triplet_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetTriplet_r_wrp
  PUBLIC :: GetTripletValues_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetTriplet_c_wrp
  PUBLIC :: GetTripletValues_c_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  PURE SUBROUTINE SetTriplet_r_wrp(ih_this,index_column,index_row,point_value) &
       & bind(c,name="SetTriplet_r_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    REAL(NTREAL), INTENT(IN) :: point_value
    TYPE(Triplet_r_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%SetTriplet(index_column, index_row, point_value)
  END SUBROUTINE SetTriplet_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  PURE SUBROUTINE GetTripletValues_r_wrp(ih_this,index_column,index_row, &
       & point_value) bind(c,name="GetTripletValues_r_wrp")
    INTEGER(kind=c_int), INTENT(IN)     :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)    :: index_column
    INTEGER(kind=c_int), INTENT(OUT)    :: index_row
    REAL(NTREAL), INTENT(OUT) :: point_value
    TYPE(Triplet_r_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%GetTriplet(index_column, index_row, point_value)
  END SUBROUTINE GetTripletValues_r_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the values of a triplet.
  PURE SUBROUTINE SetTriplet_c_wrp(ih_this,index_column,index_row,point_value) &
       & bind(c,name="SetTriplet_c_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: index_column
    INTEGER(kind=c_int), INTENT(IN)    :: index_row
    COMPLEX(NTCOMPLEX), INTENT(IN) :: point_value
    TYPE(Triplet_c_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%SetTriplet(index_column, index_row, point_value)
  END SUBROUTINE SetTriplet_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the values of a triplet.
  PURE SUBROUTINE GetTripletValues_c_wrp(ih_this,index_column,index_row, &
       & point_value) bind(c,name="GetTripletValues_c_wrp")
    INTEGER(kind=c_int), INTENT(IN)     :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)    :: index_column
    INTEGER(kind=c_int), INTENT(OUT)    :: index_row
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: point_value
    TYPE(Triplet_c_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL h_this%data%GetTriplet(index_column, index_row, point_value)
  END SUBROUTINE GetTripletValues_c_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletModule_wrp
