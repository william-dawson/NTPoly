!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A wrapper for the iterative solver parameters.
MODULE SolverParametersModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE SolverParametersModule
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_bool, c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the sparse matrix data type.
  TYPE, PUBLIC :: SolverParameters_wrp
     TYPE(SolverParameters_t), POINTER :: DATA
  END TYPE SolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructSolverParameters_wrp
  PUBLIC :: DestructSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetParametersConvergeDiff_wrp
  PUBLIC :: SetParametersMaxIterations_wrp
  PUBLIC :: SetParametersThreshold_wrp
  PUBLIC :: SetParametersBeVerbose_wrp
  PUBLIC :: SetParametersLoadBalance_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the iterat solver parameters.
  SUBROUTINE ConstructSolverParameters_wrp(ih_this) &
       & BIND(c,name="ConstructSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(SolverParameters_wrp) :: h_this

    ALLOCATE(h_this%DATA)
    h_this%DATA = SolverParameters_t()
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a iterative solver parameter type.
  SUBROUTINE DestructSolverParameters_wrp(ih_this) &
       & BIND(c,name="DestructSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructSolverParameters(h_this%DATA)
    DEALLOCATE(h_this%DATA)
    !ih_this = 0
  END SUBROUTINE DestructSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the convergence difference.
  SUBROUTINE SetParametersConvergeDiff_wrp(ih_this,new_value) &
       & BIND(c,name="SetParametersConvergeDiff_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersConvergeDiff(h_this%DATA,new_value)
  END SUBROUTINE SetParametersConvergeDiff_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the max iterations.
  SUBROUTINE SetParametersMaxIterations_wrp(ih_this,new_value) &
       & BIND(c,name="SetParametersMaxIterations_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersMaxIterations(h_this%DATA,new_value)
  END SUBROUTINE SetParametersMaxIterations_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  SUBROUTINE SetParametersThreshold_wrp(ih_this,new_value) &
       & BIND(c,name="SetParametersThreshold_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersThreshold(h_this%DATA,new_value)
  END SUBROUTINE SetParametersThreshold_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  SUBROUTINE SetParametersBeVerbose_wrp(ih_this,new_value) &
       & BIND(c,name="SetParametersBeVerbose_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersBeVerbose(h_this%DATA,LOGICAL(new_value))
  END SUBROUTINE SetParametersBeVerbose_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balancing permutation.
  SUBROUTINE SetParametersLoadBalance_wrp(ih_this,ih_new_value) &
       & BIND(c,name="SetParametersLoadBalance_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_new_value(SIZE_wrp)
    TYPE(SolverParameters_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_new_value

    h_this = TRANSFER(ih_this,h_this)
    h_new_value = TRANSFER(ih_new_value,h_new_value)
    CALL SetParametersLoadBalance(h_this%DATA,h_new_value%DATA)
  END SUBROUTINE SetParametersLoadBalance_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SolverParametersModule_wrp
