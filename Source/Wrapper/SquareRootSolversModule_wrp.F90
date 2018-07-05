!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the square root solvers module for calling from other languages.
MODULE SquareRootSolversModule_wrp
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE SquareRootSolversModule, ONLY : SquareRoot, InverseSquareRoot
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SquareRoot_wrp
  PUBLIC :: InverseSquareRoot_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse square root of a matrix.
  SUBROUTINE InverseSquareRoot_wrp(ih_Mat1, ih_InverseSquareRootMat, &
       & ih_solver_parameters) bind(c,name="InverseSquareRoot_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_InverseSquareRootMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_InverseSquareRootMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_InverseSquareRootMat = TRANSFER(ih_InverseSquareRootMat, &
         & h_InverseSquareRootMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL InverseSquareRoot(h_Mat1%data, h_InverseSquareRootMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE InverseSquareRoot_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root of a matrix.
  SUBROUTINE SquareRoot_wrp(ih_Mat1, ih_SquareRootMat, ih_solver_parameters) &
       & bind(c,name="SquareRoot_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Mat1(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_SquareRootMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_Mat1
    TYPE(Matrix_ps_wrp) :: h_SquareRootMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_Mat1 = TRANSFER(ih_Mat1,h_Mat1)
    h_SquareRootMat = TRANSFER(ih_SquareRootMat, h_SquareRootMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL SquareRoot(h_Mat1%data, h_SquareRootMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE SquareRoot_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SquareRootSolversModule_wrp
