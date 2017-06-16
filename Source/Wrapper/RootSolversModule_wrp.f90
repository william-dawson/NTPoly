!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the overlap solvers module for calling from other languages.
MODULE RootSolversModule_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE RootSolversModule, ONLY : ComputeRoot, ComputeInverseRoot
  USE DistributedBlockedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeRoot_wrp
  PUBLIC :: ComputeInverseRoot_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general matrix root.
  !! @param[in]  ih_InputMat the input matrix.
  !! @param[out] ih_OutputMat = InputMat^-1/root.
  !! @param[in]  root which root to compute.
  !! @param[in]  ih_solver_parameters parameters for the solver.
  SUBROUTINE ComputeRoot_wrp(ih_InputMat, ih_OutputMat, &
       & root, ih_solver_parameters) bind(c,name="ComputeRoot_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: root
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat,h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeRoot(h_InputMat%data, h_OutputMat%data, root, &
         & h_solver_parameters%data)
  END SUBROUTINE ComputeRoot_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general inverse matrix root.
  !! @param[in]  ih_InputMat the input matrix.
  !! @param[out] ih_OutputMat = InputMat^-1/root.
  !! @param[in]  root which root to compute.
  !! @param[in]  ih_solver_parameters parameters for the solver.
  SUBROUTINE ComputeInverseRoot_wrp(ih_InputMat, ih_OutputMat, &
       & root, ih_solver_parameters) bind(c,name="ComputeInverseRoot_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: root
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat,h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeInverseRoot(h_InputMat%data, h_OutputMat%data, root, &
         & h_solver_parameters%data)
  END SUBROUTINE ComputeInverseRoot_wrp
END MODULE RootSolversModule_wrp
