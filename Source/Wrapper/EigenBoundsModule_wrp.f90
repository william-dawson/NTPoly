!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the matrix inversion module for calling from other languages.
MODULE EigenBoundsModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE EigenBoundsModule, ONLY : GershgorinBounds, PowerBounds
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds_wrp
  PUBLIC :: PowerBounds_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !! @param[in]  ih_this Matrix.
  !! @param[out] min_value min bounds on eigenspectrum.
  !! @param[out] max_value max bounds on eigenspectrum.
  SUBROUTINE GershgorinBounds_wrp(ih_this, min_value,max_value) &
       & bind(c,name="GershgorinBounds_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(inout) :: min_value
    REAL(NTREAL), INTENT(inout) :: max_value
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)

    CALL GershgorinBounds(h_this%data, min_value, max_value)
  END SUBROUTINE GershgorinBounds_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the largestest eigenvalue using the power method.
  !! @param[in]  ih_this Matrix.
  !! @param[out] max_value max bounds on eigenspectrum.
  !! @param[in]  ih_solver_parameters parameters for the solver
  SUBROUTINE PowerBounds_wrp(ih_this, max_value, ih_solver_parameters) &
       & bind(c,name="PowerBounds_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(inout) :: max_value
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PowerBounds(h_this%data, max_value, h_solver_parameters%data)
  END SUBROUTINE PowerBounds_wrp
END MODULE EigenBoundsModule_wrp
