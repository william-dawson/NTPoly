!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A wrapper for the iterative solver parameters.
MODULE SolverParametersModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE SolverParametersModule, ONLY : SolverParameters_t, &
       & SetParametersConvergeDiff, SetParametersMaxIterations, &
       & SetParametersThreshold, SetParametersBeVerbose, &
       & SetParametersLoadBalance, SetParametersDACBaseSize, &
       & SetParametersDACBaseSparsity, DestructSolverParameters
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
  PUBLIC :: SetParametersDACBaseSize_wrp
  PUBLIC :: SetParametersDACBaseSparsity_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the iterat solver parameters.
  PURE SUBROUTINE ConstructSolverParameters_wrp(ih_this) &
       & bind(c,name="ConstructSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(SolverParameters_wrp) :: h_this

    ALLOCATE(h_this%data)
    h_this%data = SolverParameters_t()
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a iterative solver parameter type.
  PURE SUBROUTINE DestructSolverParameters_wrp(ih_this) &
       & bind(c,name="DestructSolverParameters_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructSolverParameters(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructSolverParameters_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the convergence difference.
  PURE SUBROUTINE SetParametersConvergeDiff_wrp(ih_this,new_value) &
       & bind(c,name="SetParametersConvergeDiff_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersConvergeDiff(h_this%data,new_value)
  END SUBROUTINE SetParametersConvergeDiff_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the max iterations.
  PURE SUBROUTINE SetParametersMaxIterations_wrp(ih_this,new_value) &
       & bind(c,name="SetParametersMaxIterations_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersMaxIterations(h_this%data,new_value)
  END SUBROUTINE SetParametersMaxIterations_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  PURE SUBROUTINE SetParametersThreshold_wrp(ih_this,new_value) &
       & bind(c,name="SetParametersThreshold_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersThreshold(h_this%data,new_value)
  END SUBROUTINE SetParametersThreshold_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  PURE SUBROUTINE SetParametersBeVerbose_wrp(ih_this,new_value) &
       & bind(c,name="SetParametersBeVerbose_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersBeVerbose(h_this%data,LOGICAL(new_value))
  END SUBROUTINE SetParametersBeVerbose_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balancing permutation.
  PURE SUBROUTINE SetParametersLoadBalance_wrp(ih_this,ih_new_value) &
       & bind(c,name="SetParametersLoadBalance_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_new_value(SIZE_wrp)
    TYPE(SolverParameters_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_new_value

    h_this = TRANSFER(ih_this,h_this)
    h_new_value = TRANSFER(ih_new_value,h_new_value)
    CALL SetParametersLoadBalance(h_this%data,h_new_value%data)
  END SUBROUTINE SetParametersLoadBalance_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the divide and conquer base size.
  PURE SUBROUTINE SetParametersDACBaseSize_wrp(ih_this,new_value) &
       & bind(c,name="SetParametersDACBaseSize_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersDACBaseSize(h_this%data,new_value)
  END SUBROUTINE SetParametersDACBaseSize_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the divide and conquer base sparsity.
  PURE SUBROUTINE SetParametersDACBaseSparsity_wrp(ih_this,new_value) &
       & bind(c,name="SetParametersDACBaseSparsity_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(kind=NTREAL), INTENT(IN) :: new_value
    TYPE(SolverParameters_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL SetParametersDACBaseSparsity(h_this%data,new_value)
  END SUBROUTINE SetParametersDACBaseSparsity_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SolverParametersModule_wrp
