!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the eigenvalues of a matrix.
MODULE EigenSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE DMatrixModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructMatrixDFromS, &
       & ConstructMatrixSFromD, DestructMatrix
  USE EigenBoundsModule, ONLY : GershgorinBounds, PowerBounds
#if EIGENEXA
  USE EigenExaModule, ONLY : EigenExa_s
#endif
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteComment, WriteElement
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix, MatrixMultiply, &
       & ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, GatherMatrixToProcess, &
       & FillMatrixFromTripletList, ConstructEmptyMatrix, ConvertMatrixToReal, &
       & DestructMatrix, CopyMatrix, GetMatrixTripletList, TransposeMatrix, &
       & ConjugateMatrix, FillMatrixIdentity, WriteMatrixToMatrixMarket
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, MatrixToTripletList, &
       & DestructMatrix
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & ConstructTripletList, DestructTripletList
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition
  PUBLIC :: DenseMatrixFunction
  PUBLIC :: EstimateGap
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigendecomposition of a matrix.
  !! Uses a dense routine.
  SUBROUTINE EigenDecomposition(this, eigenvalues, eigenvectors_in, nvals_in, &
       & solver_parameters_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Diagonal matrix of eigenvalues.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> The eigenvectors of a matrix.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvectors_in
    !> The number of desired eigenvalues.
    INTEGER, INTENT(IN), OPTIONAL :: nvals_in
    !> Parameters for computing
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! For Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    INTEGER :: nvals

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF
    IF (PRESENT(nvals_in)) THEN
       nvals = nvals_in
    ELSE
       nvals = this%actual_matrix_dimension
    END IF

#if EIGENEXA
    IF (PRESENT(eigenvectors_in)) THEN
       CALL EigenExa_s(this, eigenvalues, nvals, eigenvectors_in, &
            & solver_parameters_in = params)
    ELSE
       CALL EigenExa_s(this, eigenvalues, nvals, &
            & solver_parameters_in = params)
    END IF
#else
    IF (PRESENT(eigenvectors_in)) THEN
       CALL EigenSerial(this, eigenvalues, nvals, params, eigenvectors_in)
    ELSE
       CALL EigenSerial(this, eigenvalues, nvals, params)
    END IF
#endif

    !! Cleanup
    CALL DestructSolverParameters(params)
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
         REAL(KIND = NTREAL), INTENT(IN) :: val
         !> The transformed value.
         REAL(KIND = NTREAL) :: outval
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
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    !! Perform the eigendecomposition
    CALL EigenDecomposition(this, vals, solver_parameters_in = params, &
         & eigenvectors_in = vecs)

    !! Convert to a triplet list, map the triplet list, fill.
    CALL GetMatrixTripletList(vals, tlist)
    DO II = 1, tlist%CurrentSize
       tlist%DATA(II)%point_value = func(tlist%DATA(II)%point_value)
    END DO

    !! Fill
    CALL ConstructEmptyMatrix(ResultMat, this)
    CALL FillMatrixFromTripletList(ResultMat, tlist, preduplicated_in = .TRUE.)

    !! Multiply Back Together
    CALL MatrixMultiply(vecs, ResultMat, temp, threshold_in = params%threshold)
    CALL TransposeMatrix(vecs, vecsT)
    CALL ConjugateMatrix(vecsT)
    CALL MatrixMultiply(temp, vecsT, ResultMat, threshold_in = params%threshold)

    !! Cleanup
    CALL DestructMatrix(vecs)
    CALL DestructMatrix(temp)
    CALL DestructMatrix(vals)
    CALL DestructTripletList(tlist)
    CALL DestructSolverParameters(params)

  END SUBROUTINE DenseMatrixFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Estimate the HOMO-LUMO gap of a matrix.
  SUBROUTINE EstimateGap(H, K, chemical_potential, gap,  solver_parameters_in)
    !> The matrix to compute the gap of.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The density matrix.
    TYPE(Matrix_ps), INTENT(IN) :: K
    !> The chemical potential (estimated from purification)
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    !> The computed estimate of the HOMO-LUMO gap
    REAL(NTREAL), INTENT(OUT) :: gap
    !> Parameters for computing
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! For Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: KH, ShiftH
    TYPE(MatrixMemoryPool_p) :: pool
    REAL(NTREAL) :: e_min, e_max

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Estimate Gap")
       CALL EnterSubLog
       CALL PrintParameters(params)
    END IF

    !! Use the power method to get the smallest eigenvalue. This works
    !! because we project out the rydberg states, but if not we fall
    !! back to the Gershgorin Bounds
    CALL MatrixMultiply(K, H, KH, &
         & threshold_in = params%threshold, memory_pool_in = pool)
    CALL WriteElement("Estimate Minimum")
    CALL EnterSubLog
    CALL PowerBounds(KH, e_min, params)
    CALL ExitSubLog
    IF (e_min .GT. 0_NTREAL) THEN
       CALL WriteComment("Switching to GershgorinBounds")
       CALL GershgorinBounds(H, e_min, e_max)
    END IF
    CALL WriteElement("Estimated e_min", VALUE = e_min)

    !! Shift the matrix so the HOMO is the largest magnitude eigenvalue
    CALL ConstructEmptyMatrix(ShiftH, H)
    CALL FillMatrixIdentity(ShiftH)
    CALL ScaleMatrix(ShiftH, -e_min)
    CALL IncrementMatrix(H, ShiftH)

    !! Project Out the Virtual Orbitals
    CALL MatrixMultiply(K, ShiftH, KH, &
         & threshold_in = params%threshold, memory_pool_in = pool)

    !! Use the power method to get the largest eigenvalue
    CALL PowerBounds(KH, e_max, params)
    e_max = e_max + e_min
    CALL WriteElement("HOMO Estimate", VALUE = e_max)
    gap = 2.0_NTREAL * (chemical_potential - e_max)
    CALL WriteElement("Gap Estimate", VALUE = gap)

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructMatrix(KH)
    CALL DestructMatrix(ShiftH)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)

  END SUBROUTINE EstimateGap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve.
  SUBROUTINE EigenSerial(this, eigenvalues, nvals, solver_params, &
       & eigenvectors_in)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> The number of vals to compute.
    INTEGER, INTENT(IN) :: nvals
    !> The solve parameters.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_params
    !> The eigenvectors of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvectors_in

    !! Write info about the solver
    IF (solver_params%be_verbose) THEN
       CALL WriteHeader("Eigen Solver")
       CALL EnterSubLog
       CALL WriteElement(key = "Method", VALUE = "LAPACK")
       CALL WriteElement(key = "NVALS", VALUE = nvals)
       CALL ExitSubLog
       CALL PrintParameters(solver_params)
    END IF

    IF (this%is_complex) THEN
       IF (PRESENT(eigenvectors_in)) THEN
          CALL EigenSerial_c(this, eigenvalues, nvals, &
               & solver_params%threshold, eigenvectors_in)
       ELSE
          CALL EigenSerial_c(this, eigenvalues, nvals, &
               & solver_params%threshold)
       END IF
    ELSE
       IF (PRESENT(eigenvectors_in)) THEN
          CALL EigenSerial_r(this, eigenvalues, nvals, &
               & solver_params%threshold, eigenvectors_in)
       ELSE
          CALL EigenSerial_r(this, eigenvalues, nvals, &
               & solver_params%threshold)
       END IF
    END IF
  END SUBROUTINE EigenSerial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve (REAL).
  SUBROUTINE EigenSerial_r(this, eigenvalues, nvals, threshold, &
       & eigenvectors_in)
    USE DMatrixModule, ONLY : EigenDecomposition
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> Number of values to compute.
    INTEGER, INTENT(IN) :: nvals
    !> Threshold
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Matrix eigenvectors.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvectors_in
    !! Local variables
    TYPE(Matrix_lsr) :: local_s, V_s, W_s
    TYPE(Matrix_ldr) :: local_d, V, W
    TYPE(TripletList_r) :: V_t, W_t

#include "eigenexa_includes/EigenSerial.f90"

  END SUBROUTINE EigenSerial_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve (COMPLEX).
  SUBROUTINE EigenSerial_c(this, eigenvalues, nvals, threshold, &
       & eigenvectors_in)
    USE DMatrixModule, ONLY : EigenDecomposition
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvalues of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> Number of values to compute.
    INTEGER, INTENT(IN) :: nvals
    !> Threshold
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Matrix eigenvectors.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvectors_in
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
