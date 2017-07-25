!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the overlap solvers module for calling from other languages.
MODULE ExponentialSolversModule_wrp
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE ExponentialSolversModule, ONLY : ComputeExponential, ComputeLogarithm
  USE FixedSolversModule_wrp, ONLY : FixedSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeExponential_wrp
  PUBLIC :: ComputeLogarithm_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix.
  !! @param[in] ih_InputMat the input matrix
  !! @param[out] ih_OutputMat = exp(ih_InputMat)
  !! @param[in] ih_solver_parameters parameters for the solver.
  SUBROUTINE ComputeExponential_wrp(ih_InputMat, ih_OutputMat, &
       & ih_solver_parameters) bind(c,name="ComputeExponential_wrp")
    INTEGER(kind=c_int), INTENT(in)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(FixedSolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeExponential(h_InputMat%data, h_OutputMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE ComputeExponential_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the logarithm of a matrix.
  !! @param[in] ih_InputMat the input matrix
  !! @param[out] ih_OutputMat = log(ih_InputMat)
  !! @param[in] ih_solver_parameters parameters for the solver.
  SUBROUTINE ComputeLogarithm_wrp(ih_InputMat, ih_OutputMat, &
       & ih_solver_parameters) bind(c,name="ComputeLogarithm_wrp")
    INTEGER(kind=c_int), INTENT(in)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp) :: h_OutputMat
    TYPE(FixedSolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeLogarithm(h_InputMat%data, h_OutputMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE ComputeLogarithm_wrp
END MODULE ExponentialSolversModule_wrp
