!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A wrapper for the chemistry solver parameters.
MODULE IterativeSolversModule_wrp
  USE IterativeSolversModule, ONLY : IterativeSolverParameters, &
       & SetIterativeConvergeDiff, SetIterativeMaxIterations, &
       & SetIterativeThreshold, SetIterativeBeVerbose, SetIterativeLoadBalance
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_bool, c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the sparse matrix data type.
  TYPE, PUBLIC :: IterativeSolverParameters_wrp
     !> Actual data.
     TYPE(IterativeSolverParameters), POINTER :: data
  END TYPE IterativeSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructIterativeSolverParameters_wrp
  PUBLIC :: DestructIterativeSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetIterativeConvergeDiff_wrp
  PUBLIC :: SetIterativeMaxIterations_wrp
  PUBLIC :: SetIterativeThreshold_wrp
  PUBLIC :: SetIterativeBeVerbose_wrp
  PUBLIC :: SetIterativeLoadBalance_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the chemistry solver parameters.
  !! @param[inout] ih_this the solver parameter.
  PURE SUBROUTINE ConstructIterativeSolverParameters_wrp(ih_this) &
       & bind(c,name="ConstructIterativeSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(IterativeSolverParameters_wrp) :: h_this

    ALLOCATE(h_this%data)
    h_this%data = IterativeSolverParameters()
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructIterativeSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a chemistry solver parameter type.
  !! @param[inout] ih_this the object to destruct.
  PURE SUBROUTINE DestructIterativeSolverParameters_wrp(ih_this) &
       & bind(c,name="DestructIterativeSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructIterativeSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the convergence difference.
  !! @param[inout] ih_this handle to this.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeConvergeDiff_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeConvergeDiff_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeConvergeDiff(h_this%data,new_value)
  END SUBROUTINE SetIterativeConvergeDiff_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the max iterations.
  !! @param[inout] ih_this handle to this.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeMaxIterations_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeMaxIterations_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeMaxIterations(h_this%data,new_value)
  END SUBROUTINE SetIterativeMaxIterations_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  !! @param[inout] ih_this handle to this.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeThreshold_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeThreshold_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeThreshold(h_this%data,new_value)
  END SUBROUTINE SetIterativeThreshold_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  !! @param[inout] ih_this handle to this.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeBeVerbose_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeBeVerbose_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(in) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeBeVerbose(h_this%data,LOGICAL(new_value))
  END SUBROUTINE SetIterativeBeVerbose_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balancing permutation.
  !! @param[inout] ih_this handle to this.
  !! @param[in] ih_new_value to set it to.
  PURE SUBROUTINE SetIterativeLoadBalance_wrp(ih_this,ih_new_value) &
       & bind(c,name="SetIterativeLoadBalance_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_new_value(SIZE_wrp)
    TYPE(IterativeSolverParameters_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_new_value

    h_this = TRANSFER(ih_this,h_this)
    h_new_value = TRANSFER(ih_new_value,h_new_value)
    CALL SetIterativeLoadBalance(h_this%data,h_new_value%data)
  END SUBROUTINE SetIterativeLoadBalance_wrp
END MODULE IterativeSolversModule_wrp