!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the overlap solvers module for calling from other languages.
MODULE TrigonometrySolversModule_wrp
  USE FixedSolversModule_wrp, ONLY : FixedSolverParameters_wrp
  USE TrigonometrySolversModule, ONLY : Sine, Cosine
  USE DistributedBlockedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: Sine_wrp
  PUBLIC :: Cosine_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the sine of a matrix.
  !! @param[in]  ih_InputMat Input.
  !! @param[out] ih_OutputMat = sin(InputMat).
  !! @param[in]  ih_solver_parameters parameters for the solver
  SUBROUTINE Sine_wrp(ih_InputMat, ih_OutputMat, ih_solver_parameters) &
       & bind(c,name="Sine_wrp")
    INTEGER(kind=c_int), INTENT(in)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp)   :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp)   :: h_OutputMat
    TYPE(FixedSolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Sine(h_InputMat%data, h_OutputMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE Sine_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the cosine of a matrix.
  !! @param[in]  ih_InputMat Input.
  !! @param[out] ih_OutputMat = cosin(InputMat).
  !! @param[in]  ih_solver_parameters parameters for the solver
  SUBROUTINE Cosine_wrp(ih_InputMat, ih_OutputMat, ih_solver_parameters) &
       & bind(c,name="Cosine_wrp")
    INTEGER(kind=c_int), INTENT(in)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp)   :: h_InputMat
    TYPE(DistributedSparseMatrix_wrp)   :: h_OutputMat
    TYPE(FixedSolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Cosine(h_InputMat%data, h_OutputMat%data, &
         & h_solver_parameters%data)
  END SUBROUTINE Cosine_wrp
END MODULE TrigonometrySolversModule_wrp
