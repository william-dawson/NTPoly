!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing estimates of the bounds of a matrix's spectrum.
MODULE EigenBoundsModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DMatrixModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructMatrixDFromS, &
       & ConstructMatrixSFromD, EigenDecomposition, MultiplyMatrix, &
       & TransposeMatrix, IncrementMatrix, ConstructEmptyMatrix, MatrixNorm, &
       & CopyMatrix, ConjugateMatrix, DestructMatrix
  USE LinearSolversModule, ONLY : PivotedCholeskyDecomposition
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader
  USE MatrixMapsModule, ONLY : MapMatrix_psr
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, DotMatrix, &
       & IncrementMatrix, ScaleMatrix, MatrixTrace
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, GetMatrixTripletList, FillMatrixFromTripletList, &
       & ConvertMatrixToComplex, GatherMatrixToProcess, TransposeMatrix, &
       & ResizeMatrix, PrintMatrix, FillMatrixIdentity, ConjugateMatrix, &
       & ConvertMatrixToReal
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, MatrixToTripletList
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & AppendToTripletList, DestructTripletList, ConstructTripletList
  USE TripletModule, ONLY : Triplet_r
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds
  PUBLIC :: PowerBounds
  PUBLIC :: InteriorEigenvalues
  PUBLIC :: SubspaceIteration
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !> Uses Gershgorin's theorem.
  SUBROUTINE GershgorinBounds(this,min_value,max_value)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> A lower bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: min_value
    !> An uppder bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !! Local Data
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    !! Counters/Temporary
    INTEGER :: counter
    INTEGER :: local_column
    INTEGER :: ierr

    IF (this%is_complex) THEN
#define triplet_list triplet_list_c
#include "solver_includes/GershgorinBounds.f90"
#undef triplet_list
    ELSE
#define triplet_list triplet_list_r
#include "solver_includes/GershgorinBounds.f90"
#undef triplet_list
    END IF
  END SUBROUTINE GershgorinBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the maximum eigenvalue of a matrix.
  !> Uses The Power Method.
  SUBROUTINE PowerBounds(this,max_value,solver_parameters_in)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> An upper bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Data
    TYPE(Matrix_ps) :: vector, vector2, TempMat
    REAL(NTREAL) :: scale_value
    REAL(NTREAL) :: norm_value
    TYPE(TripletList_r) :: temp_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: outer_counter
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
       solver_parameters%max_iterations = 10
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Power Bounds Solver")
       CALL EnterSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Diagonal matrices serve as vectors.
    CALL ConstructEmptyMatrix(vector, this)
    CALL ConstructEmptyMatrix(vector2, this)

    !! Guess Vector
    temp_list = TripletList_r()
    IF (this%process_grid%global_rank .EQ. 0) THEN
       temp_triplet%index_row = 1
       temp_triplet%index_column = 1
       temp_triplet%point_value = 1.0_NTREAL
       CALL AppendToTripletList(temp_list,temp_triplet)
    END IF
    CALL FillMatrixFromTripletList(vector,temp_list)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", value=outer_counter-1)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL ExitSubLog
       END IF

       !! x = Ax
       CALL MatrixMultiply(this,vector,vector2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       !! x = x/||x||
       scale_value = 1.0/MatrixNorm(vector2)
       CALL ScaleMatrix(vector2,scale_value)

       !! Check if Converged
       CALL IncrementMatrix(vector2,vector,-1.0_NTREAL)
       norm_value = MatrixNorm(vector)

       CALL CopyMatrix(vector2,vector)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",value=outer_counter-1)
    END IF

    !! Compute The Largest Eigenvalue
    CALL DotMatrix(vector, vector, scale_value)
    CALL MatrixMultiply(this,vector,vector2, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL DotMatrix(vector, vector2, max_value)
    max_value = max_value / scale_value

    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Max_Eigen_Value",value=max_value)
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(vector)
    CALL DestructMatrix(vector2)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE PowerBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute interior eigenvalues of a matrix.
  SUBROUTINE InteriorEigenvalues(this, density, nel, nvals, eigvecs, &
       & eigenvalues_out, solver_parameters_in)
    !> The matrix to compute the interior eigenvalues of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The density matrix that splits the spectrum of this matrix.
    TYPE(Matrix_ps), INTENT(IN) :: density
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The number of values to compute. Negative if they should be below
    !! the gap, positive if above.
    INTEGER, INTENT(IN) :: nvals
    !> The computed eigenvectors.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigvecs
    !> The eigenvalues, stored in a sparse matrix.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local matrices
    TYPE(Matrix_ps) :: target_density
    TYPE(Matrix_ps) :: identity
    TYPE(Matrix_ps) :: vec, vecT
    TYPE(Matrix_ps) :: tempmat
    TYPE(Matrix_ps) :: reduced_matrix, shifted_matrix
    !! Local variables
    INTEGER :: target_dim
    REAL(NTREAL) :: min_val, max_val

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    !! Setup whether to compute the left or right spectrum.
    CALL ConstructEmptyMatrix(identity, density)
    CALL FillMatrixIdentity(identity)
    IF (nvals .GT. 0) THEN
       target_dim = this%actual_matrix_dimension - 0.5*nel
       CALL CopyMatrix(identity, target_density)
       CALL IncrementMatrix(density, target_density, alpha_in=-1.0_NTREAL, &
            & threshold_in=solver_parameters%threshold)
    ELSE
       target_dim = 0.5*nel
       CALL CopyMatrix(density, target_density)
    END IF

    !! Compute the pivoted cholesky decomposition.
    CALL PivotedCholeskyDecomposition(target_density, vec, target_dim, &
         & solver_parameters_in=solver_parameters)

    !! Construct the reduced matrix.
    CALL MatrixMultiply(this, vec, tempmat, &
         & threshold_in=solver_parameters%threshold)
    CALL TransposeMatrix(vec, vecT)
    IF (vec%is_complex) THEN
       CALL ConjugateMatrix(vecT)
    END IF
    CALL MatrixMultiply(vecT, tempmat, reduced_matrix, &
         & threshold_in=solver_parameters%threshold)
    CALL ResizeMatrix(reduced_matrix, target_dim)

    !! Shift the matrix to get the extreme absolute values.
    CALL GershgorinBounds(reduced_matrix, min_val, max_val)
    CALL ConstructEmptyMatrix(shifted_matrix, reduced_matrix)
    IF (nvals .GT. 0) THEN
       CALL FillMatrixIdentity(shifted_matrix)
       CALL ScaleMatrix(shifted_matrix, max_val)
       CALL IncrementMatrix(reduced_matrix, shifted_matrix, &
            & alpha_in=-1.0_NTREAL)
    ELSE
       CALL FillMatrixIdentity(shifted_matrix)
       CALL ScaleMatrix(shifted_matrix, -1.0_NTREAL*min_val)
       CALL IncrementMatrix(reduced_matrix, shifted_matrix)
    END IF

    !! Subspace iteration.
    CALL SubspaceIteration(shifted_matrix, tempmat, ABS(nvals), &
         & solver_parameters_in=solver_parameters)

    !! Rotate back
    CALL ResizeMatrix(tempmat, this%actual_matrix_dimension)
    CALL MatrixMultiply(vec, tempmat, eigvecs, &
         & threshold_in=solver_parameters%threshold)

    !! Compute eigenvalues.
    IF (PRESENT(eigenvalues_out)) THEN
       CALL TransposeMatrix(eigvecs, vecT)
       CALL MatrixMultiply(vecT, this, tempmat)
       CALL MatrixMultiply(tempmat, eigvecs, eigenvalues_out, &
            & threshold_in=solver_parameters%threshold)
       CALL CopyMatrix(eigenvalues_out, tempmat)
       CALL ResizeMatrix(tempmat, ABS(nvals))
       CALL MapMatrix_psr(tempmat, eigenvalues_out, DiagonalMap)
    END IF

    !! Cleanup
    CALL DestructMatrix(target_density)
    CALL DestructMatrix(identity)
    CALL DestructMatrix(vec)
    CALL DestructMatrix(vecT)
    CALL DestructMatrix(tempmat)
    CALL DestructMatrix(reduced_matrix)
    CALL DestructMatrix(shifted_matrix)
  END SUBROUTINE InteriorEigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute K largest eigenvalues with subspace iteration.
  SUBROUTINE SubspaceIteration(this, vecs, k, eigenvalues_out, &
       & solver_parameters_in)
    !> The matrix to compute the eigenvectors of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The eigenvectors, stored in a sparse matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: vecs
    !> The number of vectors to compute
    INTEGER, INTENT(IN) :: k
    !> The eigenvalues, stored in a sparse matrix.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Temporary matrices
    TYPE(Matrix_ps) :: temp_mat
    TYPE(Matrix_ps) :: vecs2
    TYPE(Matrix_ps) :: vecst
    TYPE(Matrix_ldr) :: ritz_values, ritz_values2
    TYPE(MatrixMemoryPool_p) :: pool
    !! Temporary triplets
    TYPE(TripletList_r) :: temp_list
    TYPE(Triplet_r) :: temp_triplet
    !! For the loop
    REAL(NTREAL) :: norm_value
    INTEGER :: II
    INTEGER :: ierr

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Subspace Iteration Solver")
       CALL EnterSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Initial guess
    CALL ConstructEmptyMatrix(temp_mat, this%actual_matrix_dimension, &
         & this%process_grid, .FALSE.)
    temp_list = TripletList_r()
    IF (this%process_grid%global_rank .EQ. 0) THEN
       DO II =1 , k
          temp_triplet%index_row = II
          temp_triplet%index_column = II
          temp_triplet%point_value = 1.0_NTREAL
          CALL AppendToTripletList(temp_list,temp_triplet)
       END DO
    END IF
    CALL FillMatrixFromTripletList(temp_mat, temp_list)
    IF (this%is_complex) THEN
       CALL ConvertMatrixToComplex(temp_mat, vecs)
    ELSE
       CALL CopyMatrix(temp_mat, vecs)
    END IF

    !! Iteration
    CALL ConstructEmptyMatrix(ritz_values, k, 1)
    ritz_values%data = 0
    DO II = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. II .GT. 1) THEN
          CALL WriteListElement(key="Round", value=II-1)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL ExitSubLog
       END IF

       !! x = Ax
       CALL MatrixMultiply(this, vecs, vecs2, memory_pool_in=pool)

       !! Orthogonalize
       IF (this%is_complex) THEN
          CALL OrthogonalizeVectors_r(vecs2, temp_mat, k, pool, &
               & solver_parameters)
       ELSE
          CALL OrthogonalizeVectors_r(vecs2, temp_mat, k, pool, &
               & solver_parameters)
       END IF

       !! Update to our next guess.
       IF (this%is_complex) THEN
          CALL UpdateIteration_c(this, temp_mat, vecs, k, ritz_values2, pool, &
               & solver_parameters)
       ELSE
          CALL UpdateIteration_r(this, temp_mat, vecs, k, ritz_values2, pool, &
               & solver_parameters)
       END IF

       !! Check convergence
       CALL IncrementMatrix(ritz_values2, ritz_values, alpha_in=-1.0_NTREAL)
       norm_value = MAXVAL(ABS(ritz_values%data))
       CALL MPI_Allreduce(MPI_IN_PLACE, norm_value, 1, MPINTREAL, MPI_SUM, &
            & this%process_grid%within_slice_comm, ierr)
       CALL CopyMatrix(ritz_values2, ritz_values)
       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",value=II-1)
    END IF

    !! Compute eigenvalues.
    IF (PRESENT(eigenvalues_out)) THEN
       CALL TransposeMatrix(vecs, vecst)
       IF (vecst%is_complex) THEN
          CALL ConjugateMatrix(vecst)
       END IF
       CALL MatrixMultiply(vecst, this, temp_mat, memory_pool_in=pool)
       CALL MatrixMultiply(temp_mat, vecs, eigenvalues_out, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL ResizeMatrix(eigenvalues_out, k)
       IF (eigenvalues_out%is_complex) THEN
          CALL ConvertMatrixToReal(eigenvalues_out, temp_mat)
       ELSE
          CALL CopyMatrix(eigenvalues_out, temp_mat)
       END IF
       CALL MapMatrix_psr(temp_mat, eigenvalues_out, DiagonalMap)
    END IF

    !! Cleanup
    CALL DestructMatrix(temp_mat)
    CALL DestructMatrix(vecs2)
    CALL DestructMatrix(vecst)
    CALL DestructMatrix(ritz_values)
    CALL DestructMatrix(ritz_values2)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructTripletList(temp_list)
  END SUBROUTINE SubspaceIteration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This helper routine will orthogonalize some vectors.
  SUBROUTINE OrthogonalizeVectors_r(invec, outvec, num_vecs, pool, params)
    !> The matrix of vectors to orthogonalize.
    TYPE(Matrix_ps), INTENT(IN) :: invec
    !> The orthogonalized vectors which are computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: outvec
    !> The number of vectors actually stored.
    INTEGER, INTENT(IN) :: num_vecs
    !> A memory pool for calculations.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !! Local variables
    TYPE(Matrix_lsr) :: local_sparse
    TYPE(Matrix_ldr) :: local_dense, local_v, local_s, local_vt
    TYPE(TripletList_r) :: triplet_list

    INCLUDE "solver_includes/OrthogonalizeVectors.f90"
  END SUBROUTINE OrthogonalizeVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This helper routine will orthogonalize some vectors.
  SUBROUTINE OrthogonalizeVectors_c(invec, outvec, num_vecs, pool, params)
    !> The matrix of vectors to orthogonalize.
    TYPE(Matrix_ps), INTENT(IN) :: invec
    !> The orthogonalized vectors which are computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: outvec
    !> The number of vectors actually stored.
    INTEGER, INTENT(IN) :: num_vecs
    !> A memory pool for calculations.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !! Local variables
    TYPE(Matrix_lsc) :: local_sparse
    TYPE(Matrix_ldc) :: local_dense, local_v, local_s, local_vt
    TYPE(TripletList_c) :: triplet_list

    INCLUDE "solver_includes/OrthogonalizeVectors.f90"
  END SUBROUTINE OrthogonalizeVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This helper routine will compute the next guess.
  !! To do this we diagonalize the matrix in the subspace of trial vectors,
  !! and rotate.
  SUBROUTINE UpdateIteration_r(mat, invec, outvec, num_vecs, ritz_values, pool, &
       & params)
    !> The matrix we are computing the eigenvalues of.
    TYPE(Matrix_ps), INTENT(IN) :: mat
    !> The matrix of trial vectors.
    TYPE(Matrix_ps), INTENT(IN) :: invec
    !> The updated guess.
    TYPE(Matrix_ps), INTENT(INOUT) :: outvec
    !> The number of vectors actually stored.
    INTEGER, INTENT(IN) :: num_vecs
    !> The ritz values are stored here
    TYPE(Matrix_ldr) :: ritz_values
    !> A memory pool for calculations.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !! Local Variables
    TYPE(Matrix_lsr) :: local_sparse
    TYPE(Matrix_ldr) :: local_dense, local_v
    TYPE(TripletList_r) :: triplet_list

    INCLUDE "solver_includes/UpdateIteration.f90"
  END SUBROUTINE UpdateIteration_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This helper routine will compute the next guess.
  !! To do this we diagonalize the matrix in the subspace of trial vectors,
  !! and rotate.
  SUBROUTINE UpdateIteration_c(mat, invec, outvec, num_vecs, ritz_values, pool, &
       & params)
    !> The matrix we are computing the eigenvalues of.
    TYPE(Matrix_ps), INTENT(IN) :: mat
    !> The matrix of trial vectors.
    TYPE(Matrix_ps), INTENT(IN) :: invec
    !> The updated guess.
    TYPE(Matrix_ps), INTENT(INOUT) :: outvec
    !> The number of vectors actually stored.
    INTEGER, INTENT(IN) :: num_vecs
    !> The ritz values are stored here
    TYPE(Matrix_ldr) :: ritz_values
    !> A memory pool for calculations.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !! Local Variables
    TYPE(Matrix_lsc) :: local_sparse
    TYPE(Matrix_ldc) :: local_dense, local_v
    TYPE(TripletList_c) :: triplet_list

    INCLUDE "solver_includes/UpdateIteration.f90"
  END SUBROUTINE UpdateIteration_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION DiagonalMap(row, column, val) RESULT(valid)
    !> The row value of an element.
    INTEGER, INTENT(INOUT) :: row
    !> The column value of an element.
    INTEGER, INTENT(INOUT) :: column
    !> The actual value of an element.
    REAL(KIND=NTREAL), INTENT(INOUT) :: val
    !> Set this to false to filter an element.
    LOGICAL :: valid

    IF (row .EQ. column) THEN
       valid = .TRUE.
    ELSE
       valid = .FALSE.
    END IF
  END FUNCTION DiagonalMap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenBoundsModule
