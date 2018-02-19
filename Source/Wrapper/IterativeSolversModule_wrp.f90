!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A wrapper for the iterative solver parameters.
MODULE IterativeSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE IterativeSolversModule, ONLY : IterativeSolverParameters_t, &
       & SetIterativeConvergeDiff, SetIterativeMaxIterations, &
       & SetIterativeThreshold, SetIterativeBeVerbose, SetIterativeLoadBalance
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_bool, c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the sparse matrix data type.
  TYPE, PUBLIC :: IterativeSolverParameters_wrp
     TYPE(IterativeSolverParameters_t), POINTER :: data
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
  !> Construct the iterat solver parameters.
  PURE SUBROUTINE ConstructIterativeSolverParameters_wrp(ih_this) &
       & bind(c,name="ConstructIterativeSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(IterativeSolverParameters_wrp) :: h_this

    ALLOCATE(h_this%data)
    h_this%data = IterativeSolverParameters_t()
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructIterativeSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a iterative solver parameter type.
  PURE SUBROUTINE DestructIterativeSolverParameters_wrp(ih_this) &
       & bind(c,name="DestructIterativeSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructIterativeSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the convergence difference.
  PURE SUBROUTINE SetIterativeConvergeDiff_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeConvergeDiff_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeConvergeDiff(h_this%data,new_value)
  END SUBROUTINE SetIterativeConvergeDiff_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the max iterations.
  PURE SUBROUTINE SetIterativeMaxIterations_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeMaxIterations_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeMaxIterations(h_this%data,new_value)
  END SUBROUTINE SetIterativeMaxIterations_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  PURE SUBROUTINE SetIterativeThreshold_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeThreshold_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeThreshold(h_this%data,new_value)
  END SUBROUTINE SetIterativeThreshold_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  PURE SUBROUTINE SetIterativeBeVerbose_wrp(ih_this,new_value) &
       & bind(c,name="SetIterativeBeVerbose_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(IN) :: new_value
    TYPE(IterativeSolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetIterativeBeVerbose(h_this%data,LOGICAL(new_value))
  END SUBROUTINE SetIterativeBeVerbose_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balancing permutation.
  PURE SUBROUTINE SetIterativeLoadBalance_wrp(ih_this,ih_new_value) &
       & bind(c,name="SetIterativeLoadBalance_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_new_value(SIZE_wrp)
    TYPE(IterativeSolverParameters_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_new_value

    h_this = TRANSFER(ih_this,h_this)
    h_new_value = TRANSFER(ih_new_value,h_new_value)
    CALL SetIterativeLoadBalance(h_this%data,h_new_value%data)
  END SUBROUTINE SetIterativeLoadBalance_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE IterativeSolversModule_wrp
