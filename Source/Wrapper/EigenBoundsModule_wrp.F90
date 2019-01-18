!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigenbounds module for calling from other languages.
MODULE EigenBoundsModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds, PowerBounds, SubspaceIteration
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds_wrp
  PUBLIC :: PowerBounds_wrp
  PUBLIC :: SubspaceIteration_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  SUBROUTINE GershgorinBounds_wrp(ih_this, min_value,max_value) &
       & bind(c,name="GershgorinBounds_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(INOUT) :: min_value
    REAL(NTREAL), INTENT(INOUT) :: max_value
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)

    CALL GershgorinBounds(h_this%data, min_value, max_value)
  END SUBROUTINE GershgorinBounds_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the largest eigenvalue using the power method.
  SUBROUTINE PowerBounds_wrp(ih_this, max_value, ih_solver_parameters) &
       & bind(c,name="PowerBounds_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(INOUT) :: max_value
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PowerBounds(h_this%data, max_value, h_solver_parameters%data)
  END SUBROUTINE PowerBounds_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute K largest eigenvalues with subspace iteration.
  SUBROUTINE SubspaceIteration_wrp(ih_this, ih_vecs, k, ih_solver_parameters) &
       & bind(c,name="SubspaceIteration_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_vecs(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: k
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_vecs
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_vecs = TRANSFER(ih_vecs,h_vecs)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SubspaceIteration(h_this%data, h_vecs%data, k, &
         & h_solver_parameters%data)
  END SUBROUTINE SubspaceIteration_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenBoundsModule_wrp
