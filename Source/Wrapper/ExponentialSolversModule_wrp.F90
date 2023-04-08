!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the exponential solvers module for calling from other languages.
MODULE ExponentialSolversModule_wrp
  USE ExponentialSolversModule
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE SolverParametersModule_wrp, ONLY : SolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeExponential_wrp
  PUBLIC :: ComputeDenseExponential_wrp
  PUBLIC :: ComputeLogarithm_wrp
  PUBLIC :: ComputeDenseLogarithm_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix.
  SUBROUTINE ComputeExponential_wrp(ih_InputMat, ih_OutputMat, &
       & ih_solver_parameters) BIND(c,name="ComputeExponential_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeExponential(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ComputeExponential_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix (dense version).
  SUBROUTINE ComputeDenseExponential_wrp(ih_InputMat, ih_OutputMat, &
       & ih_solver_parameters) BIND(c,name="ComputeDenseExponential_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeDenseExponential(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ComputeDenseExponential_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the exponential of a matrix with the pade approximation.
  SUBROUTINE ComputeExponentialPade_wrp(ih_InputMat, ih_OutputMat, &
       & ih_solver_parameters) BIND(c,name="ComputeExponentialPade_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeExponentialPade(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ComputeExponentialPade_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the logarithm of a matrix.
  SUBROUTINE ComputeLogarithm_wrp(ih_InputMat, ih_OutputMat, &
       & ih_solver_parameters) BIND(c,name="ComputeLogarithm_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeLogarithm(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ComputeLogarithm_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the logarithm of a matrix (dense version).
  SUBROUTINE ComputeDenseLogarithm_wrp(ih_InputMat, ih_OutputMat, &
       & ih_solver_parameters) BIND(c,name="ComputeDenseLogarithm_wrp")
    INTEGER(kind=c_int), INTENT(IN)    :: ih_InputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_OutputMat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN)    :: ih_solver_parameters(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_InputMat
    TYPE(Matrix_ps_wrp) :: h_OutputMat
    TYPE(SolverParameters_wrp) :: h_solver_parameters

    h_InputMat = TRANSFER(ih_InputMat,h_InputMat)
    h_OutputMat = TRANSFER(ih_OutputMat, h_OutputMat)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ComputeDenseLogarithm(h_InputMat%DATA, h_OutputMat%DATA, &
         & h_solver_parameters%DATA)
  END SUBROUTINE ComputeDenseLogarithm_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ExponentialSolversModule_wrp
