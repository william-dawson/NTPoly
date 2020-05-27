!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigenbounds module for calling from other languages.
MODULE EigenBoundsModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds, PowerBounds
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds_wrp
  PUBLIC :: PowerBounds_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  SUBROUTINE GershgorinBounds_wrp(ih_this, min_value,max_value) &
       & BIND(c,name="GershgorinBounds_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(INOUT) :: min_value
    REAL(NTREAL), INTENT(INOUT) :: max_value
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)

    CALL GershgorinBounds(h_this%DATA, min_value, max_value)
  END SUBROUTINE GershgorinBounds_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the largest eigenvalue using the power method.
  SUBROUTINE PowerBounds_wrp(ih_this, max_value, ih_solver_parameters) &
       & BIND(c,name="PowerBounds_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(INOUT) :: max_value
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PowerBounds(h_this%DATA, max_value, h_solver_parameters%DATA)
  END SUBROUTINE PowerBounds_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenBoundsModule_wrp
