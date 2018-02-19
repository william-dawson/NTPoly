!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A wrapper for the fixed solver parameters.
MODULE FixedSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE FixedSolversModule, ONLY : FixedSolverParameters_t, &
       & SetFixedThreshold, SetFixedBeVerbose, SetFixedLoadBalance
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_bool, c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the fixed solver parameters data type.
  TYPE, PUBLIC :: FixedSolverParameters_wrp
     TYPE(FixedSolverParameters_t), POINTER :: data
  END TYPE FixedSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructFixedSolverParameters_wrp
  PUBLIC :: DestructFixedSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetFixedThreshold_wrp
  PUBLIC :: SetFixedBeVerbose_wrp
  PUBLIC :: SetFixedLoadBalance_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the fixed solver parameters.
  PURE SUBROUTINE ConstructFixedSolverParameters_wrp(ih_this) &
       & bind(c,name="ConstructFixedSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(FixedSolverParameters_wrp) :: h_this

    ALLOCATE(h_this%data)
    h_this%data = FixedSolverParameters_t()
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructFixedSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a fixed solver parameter type.
  PURE SUBROUTINE DestructFixedSolverParameters_wrp(ih_this) &
       & bind(c,name="DestructFixedSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(FixedSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructFixedSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  PURE SUBROUTINE SetFixedThreshold_wrp(ih_this,new_value) &
       & bind(c,name="SetFixedThreshold_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: new_value
    TYPE(FixedSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetFixedThreshold(h_this%data,new_value)
  END SUBROUTINE SetFixedThreshold_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  PURE SUBROUTINE SetFixedBeVerbose_wrp(ih_this,new_value) &
       & bind(c,name="SetFixedBeVerbose_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(in) :: new_value
    TYPE(FixedSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetFixedBeVerbose(h_this%data,LOGICAL(new_value))
  END SUBROUTINE SetFixedBeVerbose_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balancing permutation.
  PURE SUBROUTINE SetFixedLoadBalance_wrp(ih_this,ih_new_value) &
       & bind(c,name="SetFixedLoadBalance_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_new_value(SIZE_wrp)
    TYPE(FixedSolverParameters_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_new_value

    h_this = TRANSFER(ih_this,h_this)
    h_new_value = TRANSFER(ih_new_value,h_new_value)
    CALL SetFixedLoadBalance(h_this%data,h_new_value%data)
  END SUBROUTINE SetFixedLoadBalance_wrp
END MODULE FixedSolversModule_wrp
