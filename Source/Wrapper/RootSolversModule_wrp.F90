!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the root solvers module for calling from other languages.
MODULE RootSolversModule_wrp
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE RootSolversModule, ONLY : ComputeRoot, ComputeInverseRoot
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeRoot_wrp
  PUBLIC :: ComputeInverseRoot_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general matrix root.
  SUBROUTINE ComputeRoot_wrp(ih_InputMat, ih_OutputMat, &
       & root, ih_solver_parameters) BIND(c,name="ComputeRoot_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: root
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat,h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeRoot(h_InputMat%DATA, h_OutputMat%DATA, root, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ComputeRoot_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general inverse matrix root.
  SUBROUTINE ComputeInverseRoot_wrp(ih_InputMat, ih_OutputMat, &
       & root, ih_solver_parameters) BIND(c,name="ComputeInverseRoot_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: root
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat,h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeInverseRoot(h_InputMat%DATA, h_OutputMat%DATA, root, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ComputeInverseRoot_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE RootSolversModule_wrp
