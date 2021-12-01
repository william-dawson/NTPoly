!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the eigenvalues or singular values of a matrix.
MODULE EigenSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE DMatrixModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructMatrixDFromS, &
       & ConstructMatrixSFromD, DestructMatrix
#if EIGENEXA
  USE EigenExaModule, ONLY : EigenExa_s
#endif
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, WriteElement
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply
  USE PSMatrixModule, ONLY : Matrix_ps, GatherMatrixToProcess, &
       & FillMatrixFromTripletList, ConstructEmptyMatrix, ConvertMatrixToReal, &
       & DestructMatrix, CopyMatrix, GetMatrixTripletList, TransposeMatrix, &
       & ConjugateMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, MatrixToTripletList, &
       & DestructMatrix
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & ConstructTripletList, DestructTripletList
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition
  PUBLIC :: DenseMatrixFunction
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigendecomposition of a matrix.
  !! Uses a dense routine.
  SUBROUTINE EigenDecomposition(this, eigenvectors, eigenvalues, &
       & solver_parameters_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvectors of a matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> Diagonal matrix of eigenvalues.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> Parameters for computing
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! For Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

#if EIGENEXA
    CALL EigenExa_s(this, eigenvectors, eigenvalues, params)
#else
    CALL EigenSerial(this, eigenvectors, eigenvalues, params)
#endif

  END SUBROUTINE EigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Apply an arbitrary matrix function defined by a matrix map as a
  !! transformation of the eigenvalues.
  SUBROUTINE DenseMatrixFunction(this, ResultMat, func, solver_parameters_in)
    !> The matrix to apply the function to.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The transformed matrix
    TYPE(Matrix_ps), INTENT(INOUT) :: ResultMat
    INTERFACE
       !> The procedure to apply to each eigenvalue.
       FUNCTION func(val) RESULT(outval)
         USE DataTypesModule, ONLY : NTREAL
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(IN) :: val
         !> The transformed value.
         REAL(KIND=NTREAL) :: outval
       END FUNCTION func
    END INTERFACE
    !> Parameters for computing
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! For Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: vecs, vecsT, vals, temp
    TYPE(TripletList_r) :: tlist
    INTEGER :: II

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    !! Perform the eigendecomposition
    CALL EigenDecomposition(this, vecs, vals, params)

    !! Convert to a triplet list, map the triplet list, fill.
    CALL GetMatrixTripletList(vals, tlist)
    DO II = 1, tlist%CurrentSize
       tlist%DATA(II)%point_value = func(tlist%DATA(II)%point_value)
    END DO

    !! Fill
    CALL ConstructEmptyMatrix(ResultMat, this)
    CALL FillMatrixFromTripletList(ResultMat, tlist, preduplicated_in=.TRUE.)

    !! Multiply Back Together
    CALL MatrixMultiply(vecs, ResultMat, temp, threshold_in=params%threshold)
    CALL TransposeMatrix(vecs, vecsT)
    CALL ConjugateMatrix(vecsT)
    CALL MatrixMultiply(temp, vecsT, ResultMat, threshold_in=params%threshold)

    !! Cleanup
    CALL DestructMatrix(vecs)
    CALL DestructMatrix(vals)
    CALL DestructTripletList(tlist)

  END SUBROUTINE DenseMatrixFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve.
  SUBROUTINE EigenSerial(this, eigenvectors, eigenvalues, solver_params)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvectors of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> The solve parameters.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_params

    !! Write info about the solver
    IF (solver_params%be_verbose) THEN
       CALL WriteHeader("Eigen Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="LAPACK")
       CALL ExitSubLog
       CALL PrintParameters(solver_params)
    END IF

    IF (this%is_complex) THEN
       CALL EigenSerial_c(this, eigenvectors, eigenvalues, &
            & solver_params%threshold)
    ELSE
       CALL EigenSerial_r(this, eigenvectors, eigenvalues, &
            & solver_params%threshold)
    END IF
  END SUBROUTINE EigenSerial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve (REAL).
  SUBROUTINE EigenSerial_r(this, eigenvectors, eigenvalues, threshold)
    USE DMatrixModule, ONLY : EigenDecomposition
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Matrix eigenvectors.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> Threshold
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local variables
    TYPE(Matrix_lsr) :: local_s, V_s, W_s
    TYPE(Matrix_ldr) :: local_d, V, W
    TYPE(TripletList_r) :: V_t, W_t

#include "eigenexa_includes/EigenSerial.f90"

  END SUBROUTINE EigenSerial_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve (COMPLEX).
  SUBROUTINE EigenSerial_c(this, eigenvectors, eigenvalues, threshold)
    USE DMatrixModule, ONLY : EigenDecomposition
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Matrix eigenvectors.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> Threshold
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local variables
    TYPE(Matrix_lsc) :: local_s, V_s, W_s
    TYPE(Matrix_ldc) :: local_d, V, W
    TYPE(TripletList_c) :: V_t, W_t
    TYPE(Matrix_ps) :: eigenvalues_r

#include "eigenexa_includes/EigenSerial.f90"

    CALL ConvertMatrixToReal(eigenvalues, eigenvalues_r)
    CALL CopyMatrix(eigenvalues_r, eigenvalues)
    CALL DestructMatrix(eigenvalues_r)

  END SUBROUTINE EigenSerial_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule
