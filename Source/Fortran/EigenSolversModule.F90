!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the eigenvalues or singular values of a matrix.
MODULE EigenSolversModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DMatrixModule, ONLY : Matrix_ldr, Matrix_ldc, DestructMatrix, &
       & ConstructMatrixSFromD, ConstructMatrixDFromS, EigenDecomposition
#if EIGENEXA
  USE EigenExaModule, ONLY : EigenExa_s
#endif
  USE LinearSolversModule, ONLY : PivotedCholeskyDecomposition
  USE PermutationModule, ONLY : Permutation_t, ConstructRandomPermutation, &
       & DestructPermutation
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, DestructMatrix, &
       & CopyMatrix, ConjugateMatrix, TransposeMatrix, GetMatrixSize, &
       & FillMatrixFromTripletList, CommSplitMatrix, ConvertMatrixToReal, &
       & FillMatrixIdentity, GetMatrixTripletList, PrintMatrixInformation, &
       & GatherMatrixToProcess
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, MatrixToTripletList, &
       & ConstructMatrixFromTripletList
  USE SignSolversModule, ONLY : PolarDecomposition
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & AppendToTripletList, ConstructTripletList, DestructTripletList, &
       & GetTripletAt, RedistributeTripletLists, SortTripletList
  USE TripletModule, ONLY : Triplet_r
  USE MPI
  IMPLICIT NONE

  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReferenceEigenDecomposition
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine uses a dense eigensolver for reference purposes.
  SUBROUTINE ReferenceEigenDecomposition(this, eigenvectors, eigenvalues_in, &
       & solver_parameters_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvectors of a matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> Diagonal matrix of eigenvalues.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_in
    !> Parameters for computing
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! For Handling Optional Parameters
    TYPE(SolverParameters_t) :: fixed_params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       fixed_params = solver_parameters_in
    ELSE
       fixed_params = SolverParameters_t()
    END IF

#if EIGENEXA
    IF (this%is_complex) THEN
       IF (PRESENT(eigenvalues_in)) THEN
          CALL EigenSerial(this, eigenvectors, eigenvalues_out=eigenvalues_in, &
               & solver_parameters_in=fixed_params)
       ELSE
          CALL EigenSerial(this, eigenvectors, solver_parameters_in=fixed_params)
       END IF
    ELSE
       IF (PRESENT(eigenvalues_in)) THEN
          CALL EigenExa_s(this, eigenvectors, eigenvalues_out=eigenvalues_in, &
               & solver_parameters_in=fixed_params)
       ELSE
          CALL EigenExa_s(this, eigenvectors, solver_parameters_in=fixed_params)
       ENDIF
    END IF
#else
    IF (PRESENT(eigenvalues_in)) THEN
       CALL EigenSerial(this, eigenvectors, eigenvalues_out=eigenvalues_in, &
            & solver_parameters_in=fixed_params)
    ELSE
       CALL EigenSerial(this, eigenvectors, solver_parameters_in=fixed_params)
    END IF
#endif

  END SUBROUTINE ReferenceEigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve.
  SUBROUTINE EigenSerial(this, eigenvectors, eigenvalues_out, &
       & solver_parameters_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvectors of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
    !> The solve parameters.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Local Data
    TYPE(SolverParameters_t) :: fixed_params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       fixed_params = solver_parameters_in
    ELSE
       fixed_params = SolverParameters_t()
    END IF

    IF (this%is_complex) THEN
       IF (PRESENT(eigenvalues_out)) THEN
          CALL EigenSerial_c(this, eigenvectors, fixed_params, eigenvalues_out)
       ELSE
          CALL EigenSerial_c(this, eigenvectors, fixed_params)
       END IF
    ELSE
       IF (PRESENT(eigenvalues_out)) THEN
          CALL EigenSerial_r(this, eigenvectors, fixed_params, eigenvalues_out)
       ELSE
          CALL EigenSerial_r(this, eigenvectors, fixed_params)
       END IF
    END IF
  END SUBROUTINE EigenSerial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve (REAL).
  SUBROUTINE EigenSerial_r(this, eigenvectors, fixed_params, eigenvalues_out)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Matrix eigenvectors.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> The solve parameters.
    TYPE(SolverParameters_t), INTENT(IN) :: fixed_params
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
    !! Local Data
    TYPE(TripletList_r) :: triplet_w, triplet_list
    TYPE(Matrix_lsr) :: sparse
    TYPE(Matrix_lsr) :: local_a, local_v
    TYPE(Matrix_ldr) :: dense_a, dense_v, dense_w

    INCLUDE "solver_includes/EigenSerial.f90"
  END SUBROUTINE EigenSerial_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve (COMPLEX).
  SUBROUTINE EigenSerial_c(this, eigenvectors, fixed_params, eigenvalues_out)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Matrix eigenvectors.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> The solve parameters.
    TYPE(SolverParameters_t), INTENT(IN) :: fixed_params
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
    !! Local Data
    TYPE(TripletList_c) :: triplet_w, triplet_list
    TYPE(Matrix_lsc) :: sparse
    TYPE(Matrix_lsc) :: local_a, local_v
    TYPE(Matrix_ldc) :: dense_a, dense_v, dense_w

    INCLUDE "solver_includes/EigenSerial.f90"
  END SUBROUTINE EigenSerial_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule
