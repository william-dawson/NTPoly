!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A wrapper for the process grid module.
MODULE ProcessGridModule_wrp
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid, &
       & GetMySlice, GetMyColumn, GetMyRow, global_grid
  USE ISO_C_BINDING, ONLY : c_int, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructProcessGrid_wrp
  PUBLIC :: DestructProcessGrid_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the process grid construction routine.
  SUBROUTINE ConstructProcessGrid_wrp(world_comm_, process_rows_, &
       & process_columns_, process_slices_, be_verbose) &
       & bind(c,name="ConstructProcessGrid_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: world_comm_
    INTEGER(kind=c_int), INTENT(IN) :: process_rows_
    INTEGER(kind=c_int), INTENT(IN) :: process_columns_
    INTEGER(kind=c_int), INTENT(IN) :: process_slices_
    LOGICAL(kind=c_bool), INTENT(IN) :: be_verbose
    CALL ConstructProcessGrid(world_comm_, process_rows_, process_columns_, &
         & process_slices_, be_verbose)
  END SUBROUTINE ConstructProcessGrid_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the slice of the current process.
  FUNCTION GetMySlice_wrp() RESULT(return_val) bind(c,name="GetMySlice_wrp")
    !! Parameters
    INTEGER(kind=c_int) :: return_val
    return_val = GetMySlice()
  END FUNCTION GetMySlice_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the column of the current process.
  FUNCTION GetMyColumn_wrp() RESULT(return_val) bind(c,name="GetMyColumn_wrp")
    !! Parameters
    INTEGER(kind=c_int) :: return_val
    return_val = GetMyColumn()
  END FUNCTION GetMyColumn_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the row of the current process.
  FUNCTION GetMyRow_wrp() RESULT(return_val) bind(c,name="GetMyRow_wrp")
    !! Parameters
    INTEGER(kind=c_int) :: return_val
    return_val = GetMyRow()
  END FUNCTION GetMyRow_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the process grid construction routine.
  SUBROUTINE DestructProcessGrid_wrp() bind(c,name="DestructProcessGrid_wrp")
    CALL DestructProcessGrid(global_grid)
  END SUBROUTINE DestructProcessGrid_wrp
END MODULE ProcessGridModule_wrp