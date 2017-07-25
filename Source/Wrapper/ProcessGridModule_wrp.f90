!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A wrapper for the process grid module.
MODULE ProcessGridModule_wrp
  USE ProcessGridModule, ONLY : ConstructProcessGrid
  USE ISO_C_BINDING, ONLY : c_int, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructProcessGrid_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the process grid construction routine.
  !> @param[in] world_comm_ a global communicator to split.
  !> @param[in] process_rows_ number of rows.
  !> @param[in] process_columns_ number of columns.
  !> @param[in] process_slices_ number of slices.
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
END MODULE ProcessGridModule_wrp
