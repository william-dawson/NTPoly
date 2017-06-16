!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the load balancer module for calling from other languages.
MODULE PermutationModule_wrp
  USE PermutationModule, ONLY : Permutation_t, ConstructDefaultPermutation, &
       ConstructReversePermutation, ConstructRandomPermutation, DestructPermutation
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the permutation data type.
  TYPE, PUBLIC :: Permutation_wrp
     !> Actual data.
     TYPE(Permutation_t), POINTER :: data
  END TYPE Permutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructDefaultPermutation_wrp
  PUBLIC :: ConstructReversePermutation_wrp
  PUBLIC :: ConstructRandomPermutation_wrp
  PUBLIC :: DestructPermutation_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that preserves the original order.
  !! @param[inout] ih_this handle to the permutation to construct.
  !! @param[in] matrix_dimension size of the matrix.
  SUBROUTINE ConstructDefaultPermutation_wrp(ih_this, matrix_dimension) &
       & bind(c,name="ConstructDefaultPermutation_wrp")
    INTEGER(kind=c_int), INTENT(out) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: matrix_dimension
    TYPE(Permutation_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructDefaultPermutation(h_this%data, matrix_dimension)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructDefaultPermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that reverses the original order.
  !! @param[inout] ih_this handle to the permutation to construct.
  !! @param[in] matrix_dimension size of the matrix.
  SUBROUTINE ConstructReversePermutation_wrp(ih_this, matrix_dimension) &
       & bind(c,name="ConstructReversePermutation_wrp")
    INTEGER(kind=c_int), INTENT(out) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: matrix_dimension
    TYPE(Permutation_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructReversePermutation(h_this%data, matrix_dimension)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructReversePermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that has a random order.
  !! @param[inout] ih_this handle to the permutation to construct.
  !! @param[in] matrix_dimension size of the matrix.
  SUBROUTINE ConstructRandomPermutation_wrp(ih_this, matrix_dimension) &
       & bind(c,name="ConstructRandomPermutation_wrp")
    INTEGER(kind=c_int), INTENT(out) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: matrix_dimension
    TYPE(Permutation_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructRandomPermutation(h_this%data, matrix_dimension)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructRandomPermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a permutation object.
  !! @param[inout] ih_this handle to the permutation to destruct.
  PURE SUBROUTINE DestructPermutation_wrp(ih_this) &
       & bind(c,name="DestructPermutation_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(Permutation_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructPermutation(h_this%data)
  END SUBROUTINE DestructPermutation_wrp
END MODULE PermutationModule_wrp
