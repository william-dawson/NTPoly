!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the permutation module for calling from other languages.
MODULE PermutationModule_wrp
  USE PermutationModule, ONLY : Permutation_t, ConstructDefaultPermutation, &
       & ConstructReversePermutation, ConstructRandomPermutation, &
       & DestructPermutation
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the permutation data type.
  TYPE, PUBLIC :: Permutation_wrp
     TYPE(Permutation_t), POINTER :: DATA
  END TYPE Permutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructDefaultPermutation_wrp
  PUBLIC :: ConstructReversePermutation_wrp
  PUBLIC :: ConstructRandomPermutation_wrp
  PUBLIC :: DestructPermutation_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that preserves the original order.
  SUBROUTINE ConstructDefaultPermutation_wrp(ih_this, matrix_dimension) &
       & BIND(c,name="ConstructDefaultPermutation_wrp")
    INTEGER(kind=c_int), INTENT(OUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: matrix_dimension
    TYPE(Permutation_wrp) :: h_this

    ALLOCATE(h_this%DATA)
    CALL ConstructDefaultPermutation(h_this%DATA, matrix_dimension)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructDefaultPermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that reverses the original order.
  SUBROUTINE ConstructReversePermutation_wrp(ih_this, matrix_dimension) &
       & BIND(c,name="ConstructReversePermutation_wrp")
    INTEGER(kind=c_int), INTENT(OUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: matrix_dimension
    TYPE(Permutation_wrp) :: h_this

    ALLOCATE(h_this%DATA)
    CALL ConstructReversePermutation(h_this%DATA, matrix_dimension)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructReversePermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that has a random order.
  SUBROUTINE ConstructRandomPermutation_wrp(ih_this, matrix_dimension) &
       & BIND(c,name="ConstructRandomPermutation_wrp")
    INTEGER(kind=c_int), INTENT(OUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: matrix_dimension
    TYPE(Permutation_wrp) :: h_this

    ALLOCATE(h_this%DATA)
    CALL ConstructRandomPermutation(h_this%DATA, matrix_dimension)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructRandomPermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a permutation object.
  SUBROUTINE DestructPermutation_wrp(ih_this) &
       & BIND(c,name="DestructPermutation_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(Permutation_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructPermutation(h_this%DATA)
  END SUBROUTINE DestructPermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PermutationModule_wrp
