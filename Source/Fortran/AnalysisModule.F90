!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Methods for analyzing the results of electronic structure calculations.
MODULE AnalysisModule
  USE CholeskyModule, ONLY : ConstructRankLookup, AppendToVector, &
       & BroadcastVector, ConstructDiag, DotAllHelper, DotAllPivoted, &
       & GatherMatrixColumn, GetPivot, UnpackCholesky
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DensityMatrixSolversModule, ONLY : TRS4
  USE DMatrixModule, ONLY : Matrix_ldr, DestructMatrix, &
       & ConstructMatrixDFromS
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteHeader, WriteListElement
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, SimilarityTransform
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & TransposeMatrix, DestructMatrix, ConjugateMatrix, CopyMatrix, &
       & FillMatrixIdentity, PrintMatrixInformation, MergeMatrixLocalBlocks, &
       & GetMatrixSlice
  USE SMatrixModule, ONLY : Matrix_lsr
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PivotedCholeskyDecomposition
  PUBLIC :: ReduceDimension
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Pivoted Cholesky Decomposition of a Hermitian Semi-Definite
  !> matrix. This is one way to generate localized orbitals.
  SUBROUTINE PivotedCholeskyDecomposition(AMat, LMat, rank_in, &
       & solver_parameters_in)
    !> The matrix A, must be hermitian, positive semi-definite.
    TYPE(Matrix_ps), INTENT(IN)  :: AMat
    !> The matrix computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: LMat
    !> The target rank of the matrix.
    INTEGER, INTENT(IN) :: rank_in
    !> Tarameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! For Pivoting
    INTEGER, DIMENSION(:), ALLOCATABLE :: pivot_vector
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: diag
    INTEGER :: pi_j
    !! Local Variables
    TYPE(Matrix_lsr) :: sparse_a
    TYPE(Matrix_lsr) :: acol
    TYPE(Matrix_ldr) :: dense_a
    !! For Storing The Results
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_column_l
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_l
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: values_l
    !! Broadcasted Column
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_index
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_values
    INTEGER :: recv_num_values
    !! Which pivots to treat locally
    INTEGER, DIMENSION(:), ALLOCATABLE :: local_pivots
    INTEGER :: num_local_pivots
    !! Root Lookups
    INTEGER, DIMENSION(:), ALLOCATABLE :: col_root_lookup
    !! Temporary Variables
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: a_buf
    INTEGER :: II, JJ, local_JJ
    INTEGER :: local_pi_i, local_pi_j
    REAL(NTREAL) :: Aval, insert_value, inverse_factor
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: dot_values
    INTEGER :: col_root
    INTEGER :: ierr

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    !! Print out parameters
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Linear Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Pivoted Cholesky Decomposition")
       CALL WriteElement(key="Target_Rank", VALUE=rank_in)
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("aquilante2006fast")
       CALL ExitSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    CALL ConstructEmptyMatrix(LMat, AMat)

    !! Construct the pivot vector
    ALLOCATE(pivot_vector(AMat%actual_matrix_dimension))
    DO JJ = 1, AMat%actual_matrix_dimension
       pivot_vector(JJ) = JJ
    END DO

    !! First get the local matrix in a dense recommendation for quick lookup
    CALL MergeMatrixLocalBlocks(AMat, sparse_a)
    CALL ConstructMatrixDFromS(sparse_a, dense_a)
    ALLOCATE(local_pivots(sparse_a%columns))

    !! Extract the diagonal
    ALLOCATE(diag(sparse_a%columns))
    CALL ConstructDiag(AMat, Lmat%process_grid, dense_a, diag)

    !! Root Lookups
    ALLOCATE(col_root_lookup(AMat%logical_matrix_dimension))
    CALL ConstructRankLookup(AMat, LMat%process_grid, col_root_lookup)

    !! Allocate space for L
    ALLOCATE(values_per_column_l(sparse_a%columns))
    ALLOCATE(index_l(sparse_a%rows, sparse_a%columns))
    ALLOCATE(values_l(sparse_a%rows, sparse_a%columns))
    values_per_column_l = 0

    !! Buffer for fast indexing of A
    ALLOCATE(a_buf(sparse_a%columns))
    a_buf = 0

    !! Allocate space for a received column
    ALLOCATE(recv_index(sparse_a%rows))
    ALLOCATE(recv_values(sparse_a%rows))
    ALLOCATE(dot_values(sparse_a%columns))

    !! Pregather the full column of A.
    CALL GatherMatrixColumn(sparse_a, acol, LMat%process_grid)

    !! Main Loop
    DO JJ = 1, rank_in
       !! Pick a pivot vector
       local_JJ  = JJ - AMat%start_row + 1
       CALL GetPivot(AMat, LMat%process_grid, JJ, pivot_vector, diag, pi_j, &
            & insert_value, local_pivots, num_local_pivots)

       !! l[pi[j],j] = sqrt(d[pi[j]])
       IF (pi_j .GE. AMat%start_column .AND. pi_j .LT. AMat%end_column) THEN
          local_pi_j = pi_j - AMat%start_column + 1
          insert_value = SQRT(insert_value)
          inverse_factor = 1.0_NTREAL/insert_value
          !! Insert
          IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
             CALL AppendToVector(values_per_column_l(local_pi_j), &
                  & index_l(:,local_pi_j), values_l(:, local_pi_j), &
                  & local_JJ, insert_value)
          END IF
          !! Extract column J for sending later
          recv_num_values = values_per_column_l(local_pi_j)
          recv_index(:recv_num_values) = index_l(:recv_num_values,local_pi_j)
          recv_values(:recv_num_values) = values_l(:recv_num_values, local_pi_j)
       END IF
       col_root = col_root_lookup(pi_j)

       !! Broadcast column JJ, and Inverse Factor
       CALL BroadcastVector(recv_num_values, recv_index, recv_values, &
            & col_root, LMat%process_grid%row_comm)
       CALL MPI_Bcast(inverse_factor, 1, MPINTREAL, col_root, &
            & LMat%process_grid%row_comm, ierr)

       !! Extract the row of A to a dense matrix for easy lookup
       DO II = MAX(acol%outer_index(pi_j),1), acol%outer_index(pi_j+1)
          a_buf(acol%inner_index(II)) = acol%values(II)
       END DO

       !! Compute Dot Products
       CALL DotAllPivoted(recv_num_values, recv_index, recv_values, &
            & values_per_column_l, index_l, values_l, local_pivots, &
            & num_local_pivots, dot_values, LMat%process_grid%column_comm)

       !! Loop over other columns
       DO II = 1, num_local_pivots
          !! Insert Into L
          local_pi_i = local_pivots(II)
          Aval = a_buf(local_pi_i)
          insert_value = inverse_factor * (Aval - dot_values(II))
          IF (ABS(insert_value) .GT. solver_parameters%threshold) THEN
             IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
                CALL AppendToVector(values_per_column_l(local_pi_i), &
                     & index_l(:,local_pi_i), values_l(:, local_pi_i), &
                     & local_JJ, insert_value)
             END IF
          END IF
          !! Update Diagonal Array
          diag(local_pi_i) = diag(local_pi_i) - insert_value**2
       END DO

       !! Clear up the A buffer
       DO II = MAX(acol%outer_index(pi_j),1), acol%outer_index(pi_j+1)
          a_buf(acol%inner_index(II)) = 0
       END DO

    END DO

    !! Unpack
    CALL UnpackCholesky(values_per_column_l, index_l, values_l, LMat)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(LMat)
       CALL ExitSubLog
    END IF
    DEALLOCATE(local_pivots)
    DEALLOCATE(pivot_vector)
    DEALLOCATE(diag)
    CALL DestructMatrix(dense_a)
    DEALLOCATE(values_per_column_l)
    DEALLOCATE(index_l)
    DEALLOCATE(values_l)
    DEALLOCATE(recv_index)
    DEALLOCATE(recv_values)
    DEALLOCATE(dot_values)
    DEALLOCATE(a_buf)
    DEALLOCATE(col_root_lookup)
    CALL DestructMatrix(sparse_a)
    CALL DestructMatrix(acol)
  END SUBROUTINE PivotedCholeskyDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> When we want to only compute the first n eigenvalues of a matrix, this
  !> routine will project out the higher eigenvalues.
  SUBROUTINE ReduceDimension(this, dim, ReducedMat, solver_parameters_in)
    !> The starting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The number of eigenvalues ot keep.
    INTEGER, INTENT(IN) :: dim
    !> a dimxdim matrix with the same first n eigenvalues as the first.
    TYPE(Matrix_ps), INTENT(INOUT) :: ReducedMat
    !> The solver parameters.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Local Variables - matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: PMat
    TYPE(Matrix_ps) :: PVec, PVecT
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: VAV
    !! Special parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
    END IF

    !! Identity matrix passed instead of ISQ
    CALL ConstructEmptyMatrix(Identity, this)
    CALL FillMatrixIdentity(Identity)

    !! Purify
    CALL TRS4(this, Identity, REAL(dim, KIND=NTREAL), PMat, &
         & solver_parameters_in=params)

    !! Compute Eigenvectors of the Density Matrix
    CALL PivotedCholeskyDecomposition(PMat, PVec, dim, &
         & solver_parameters_in=params)
    CALL TransposeMatrix(PVec, PVecT)
    IF (PVecT%is_complex) THEN
       CALL ConjugateMatrix(PVecT)
    END IF

    !! Rotate to the divided subspace
    CALL SimilarityTransform(this, PVecT, PVec, VAV, &
         & threshold_in=params%threshold)

    !! Extract
    CALL GetMatrixSlice(VAV, ReducedMat, 1, dim, 1, dim)

    CALL DestructMatrix(Identity)
    CALL DestructMatrix(PMat)
    CALL DestructMatrix(PVec)
    CALL DestructMatrix(PVecT)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(VAV)
    CALL DestructSolverParameters(params)
  END SUBROUTINE ReduceDimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE AnalysisModule
