!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the singular values of a matrix.
MODULE SingularValueSolversModule
  USE EigenSolversModule, ONLY : EigenDecomposition
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, WriteElement
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply
  USE PSMatrixModule, ONLY : Matrix_ps, DestructMatrix
  USE SignSolversModule, ONLY : PolarDecomposition
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SingularValueDecomposition
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the singular values and singular vectors of a matrix.
  SUBROUTINE SingularValueDecomposition(this, left_vectors, &
       & right_vectors, singularvalues, solver_parameters_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> A matrix containing the left singular vectors.
    TYPE(Matrix_ps), INTENT(INOUT) :: left_vectors
    !> A matrix containing the right singular vectors.
    TYPE(Matrix_ps), INTENT(INOUT) :: right_vectors
    !> A diagonal matrix containing the singularvalues.
    TYPE(Matrix_ps), INTENT(INOUT) :: singularvalues
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: UMat, HMat

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Singular Value Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Polar")
       CALL PrintParameters(params)
    END IF

    !! First compute the polar decomposition of the matrix.
    CALL PolarDecomposition(this, UMat, HMat, params)

    !! Compute the eigen decomposition of the hermitian matrix
    CALL EigenDecomposition(HMat, singularvalues, &
         & eigenvectors_in=right_vectors, solver_parameters_in=params)

    !! Compute the left singular vectors
    CALL MatrixMultiply(UMat, right_vectors, left_vectors, &
         & threshold_in=params%threshold)

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
    CALL DestructMatrix(UMat)
    CALL DestructMatrix(HMat)
  END SUBROUTINE SingularValueDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SingularValueSolversModule
