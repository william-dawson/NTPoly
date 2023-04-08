!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the trigonometry solvers module for calling from other languages.
MODULE TrigonometrySolversModule_wrp
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE TrigonometrySolversModule
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: Sine_wrp
  PUBLIC :: Cosine_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the sine of a matrix.
  SUBROUTINE Sine_wrp(ih_InputMat, ih_OutputMat, ih_solver_parameters) &
       & BIND(c,name="Sine_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp)  :: h_InputMat
    TYPE(Matrix_ps_wrp)  :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Sine(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE Sine_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the sine of a matrix (dense version).
  SUBROUTINE DenseSine_wrp(ih_InputMat, ih_OutputMat, ih_solver_parameters) &
       & BIND(c,name="DenseSine_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp)  :: h_InputMat
    TYPE(Matrix_ps_wrp)  :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL DenseSine(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE DenseSine_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the cosine of a matrix.
  SUBROUTINE Cosine_wrp(ih_InputMat, ih_OutputMat, ih_solver_parameters) &
       & BIND(c,name="Cosine_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp)  :: h_InputMat
    TYPE(Matrix_ps_wrp)  :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL Cosine(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE Cosine_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the cosine of a matrix (dense version).
  SUBROUTINE DenseCosine_wrp(ih_InputMat, ih_OutputMat, ih_solver_parameters) &
       & BIND(c,name="DenseCosine_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp)  :: h_InputMat
    TYPE(Matrix_ps_wrp)  :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL DenseCosine(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE DenseCosine_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TrigonometrySolversModule_wrp
