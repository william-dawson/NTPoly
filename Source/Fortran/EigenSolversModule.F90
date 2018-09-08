!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the eigenvalues or singular values of a matrix.
MODULE EigenSolversModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DensityMatrixSolversModule, ONLY : HPCP, TRS2
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
       & FillMatrixIdentity, GetMatrixTripletList, PrintMatrixInformation
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
  PUBLIC :: SingularValueDecomposition
  PUBLIC :: SplittingEigenDecomposition
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
  !> Compute the eigenvalues and eigenvectors of a matrix.
  SUBROUTINE SplittingEigenDecomposition(this, eigenvectors, eigenvalues, &
       & num_values_in, splits_in, solver_parameters_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A matrix containing the eigenvectors of a matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> A diagonal matrix containing the eigenvalues.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> The number of eigenvalues to compute.
    INTEGER, INTENT(IN), OPTIONAL :: num_values_in
    !> Where to split the spectrum.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: splits_in
    !> Parameters for the solve.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    INTEGER :: num_values
    INTEGER, DIMENSION(:), ALLOCATABLE :: splits
    !! Local
    TYPE(Matrix_ps) :: eigenvectorsT, TempMat, eigenvalues_r
    TYPE(Matrix_ps) :: ReducedMat
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: new_list
    TYPE(Triplet_r) :: temporary
    INTEGER :: counter

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF
    IF (PRESENT(num_values_in)) THEN
       num_values = num_values_in
    ELSE
       num_values = this%actual_matrix_dimension
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Eigen Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="recursive")
       CALL WriteElement(key="Target_Values", int_value_in=num_values)
       CALL PrintParameters(solver_parameters)
    END IF

    !! Rotate so that we are only computing the target number of values
    IF (num_values .LT. this%actual_matrix_dimension) THEN
       CALL ReduceDimension(this, num_values, solver_parameters, ReducedMat)
    ELSE
       CALL CopyMatrix(this, ReducedMat)
    END IF

    !! Compute the splits
    IF (PRESENT(splits_in)) THEN
       ALLOCATE(splits(SIZE(splits_in)))
       splits = splits_in
    ELSE
       CALL FillSplits(ReducedMat, splits)
    END IF

    !! Actual Solve
    CALL EigenRecursive(ReducedMat, eigenvectors, splits, solver_parameters)

    !! Compute the eigenvalues
    CALL TransposeMatrix(eigenvectors, eigenvectorsT)
    IF (eigenvectorsT%is_complex) THEN
       CALL ConjugateMatrix(eigenvectorsT)
    END IF
    CALL MatrixMultiply(eigenvectorsT, ReducedMat, TempMat, &
         & threshold_in=solver_parameters%threshold)
    CALL MatrixMultiply(TempMat, eigenvectors, eigenvalues, &
         & threshold_in=solver_parameters%threshold)

    !! Get rid of the off diagonal elements
    IF (eigenvalues%is_complex) THEN
       CALL ConvertMatrixToReal(eigenvalues, eigenvalues_r)
       CALL CopyMatrix(eigenvalues_r, eigenvalues)
    END IF
    CALL GetMatrixTripletList(eigenvalues,triplet_list)
    new_list = TripletList_r()
    DO counter=1,triplet_list%CurrentSize
       CALL GetTripletAt(triplet_list,counter,temporary)
       IF (temporary%index_row .EQ. temporary%index_column) THEN
          CALL AppendToTripletList(new_list,temporary)
       END IF
    END DO
    CALL DestructMatrix(this)
    CALL ConstructEmptyMatrix(this, eigenvectorsT)
    CALL FillMatrixFromTripletList(eigenvalues, new_list, &
         & preduplicated_in=.TRUE.)
    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(new_list)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    IF (ALLOCATED(splits)) THEN
       DEALLOCATE(splits)
    END IF

    CALL DestructMatrix(eigenvectorsT)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(eigenvalues_r)
    CALL DestructMatrix(ReducedMat)
  END SUBROUTINE SplittingEigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ps) :: UMat, HMat

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Singular Value Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Recursive")
       CALL PrintParameters(solver_parameters)
    END IF

    !! First compute the polar decomposition of the matrix.
    CALL PolarDecomposition(this, UMat, HMat, solver_parameters)

    !! Compute the eigen decomposition of the hermitian matrix
    CALL SplittingEigenDecomposition(HMat, right_vectors, singularvalues, &
         & solver_parameters_in=solver_parameters)

    !! Compute the left singular vectors
    CALL MatrixMultiply(UMat, right_vectors, left_vectors, &
         & threshold_in=solver_parameters%threshold)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructMatrix(UMat)
    CALL DestructMatrix(HMat)
  END SUBROUTINE SingularValueDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The recursive workhorse routine for the eigendecompositon.
  RECURSIVE SUBROUTINE EigenRecursive(this, eigenvectors, splits, &
       & solver_parameters)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A matrix containing the eigenvectors of a matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    !> Where to split the spectrum.
    INTEGER, DIMENSION(:), INTENT(IN) :: splits
    !> Parameters for the solvers.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
    !! Local Variables - matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: PMat, PHoleMat
    TYPE(Matrix_ps) :: PVec, PHoleVec
    TYPE(Matrix_ps) :: StackV, StackVT
    TYPE(Matrix_ps) :: VAV
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: LeftMat, RightMat
    TYPE(Matrix_ps) :: LeftVectors, RightVectors
    !! Local Variables - For Splitting
    INTEGER :: mat_dim, left_dim, right_dim, midpoint
    REAL(NTREAL) :: sparsity
    !! Special parameters
    TYPE(SolverParameters_t) :: balanced_it
    TYPE(Permutation_t) :: permutation

    mat_dim = this%actual_matrix_dimension
    CALL ConstructEmptyMatrix(Identity, this)
    CALL FillMatrixIdentity(Identity)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteListElement(key="Iteration_Size", int_value_in=mat_dim)
       CALL EnterSubLog
       CALL PrintMatrixInformation(this)
       CALL ExitSubLog
    END IF
    sparsity = REAL(GetMatrixSize(this),KIND=NTREAL) / &
         & (REAL(this%actual_matrix_dimension,KIND=NTREAL)**2)

    !! Base Case
    IF (SIZE(splits) .LE. 0 .OR. &
         & mat_dim .LE. solver_parameters%dac_base_size .OR. &
         & sparsity .GT. solver_parameters%dac_base_sparsity) THEN
       CALL ReferenceEigenDecomposition(this, eigenvectors, &
            & solver_parameters_in=solver_parameters)
    ELSE
       !! Setup
       midpoint = SIZE(splits)/2 + 1
       left_dim = splits(midpoint)
       right_dim = mat_dim - left_dim

       !! Setup special parameters for purification
       CALL ConstructRandomPermutation(permutation, &
            & this%logical_matrix_dimension)
       balanced_it = SolverParameters_t(&
            & converge_diff_in = solver_parameters%converge_diff, &
            & threshold_in = solver_parameters%threshold, &
            & max_iterations_in = solver_parameters%max_iterations, &
            & BalancePermutation_in=permutation)

       !! Purify
       CALL TRS2(this, Identity, left_dim*2, PMat, &
            & solver_parameters_in=balanced_it)
       CALL DestructPermutation(permutation)
       CALL CopyMatrix(Identity, PHoleMat)
       CALL IncrementMatrix(PMat, PHoleMat, alpha_in=REAL(-1.0,NTREAL), &
            & threshold_in=solver_parameters%threshold)

       !! Compute Eigenvectors of the Density Matrix
       CALL PivotedCholeskyDecomposition(PMat, PVec, left_dim, &
            & solver_parameters_in=solver_parameters)
       CALL DestructMatrix(PMat)
       CALL PivotedCholeskyDecomposition(PHoleMat, PHoleVec, right_dim, &
            & solver_parameters_in=solver_parameters)
       CALL DestructMatrix(PHoleMat)
       CALL ConstructEmptyMatrix(StackV, this)
       CALL StackMatrices(PVec, PHoleVec, left_dim, 0, StackV)
       CALL DestructMatrix(PVec)
       CALL DestructMatrix(PHoleVec)

       !! Rotate to the divided subspace
       CALL MatrixMultiply(this, StackV, TempMat, &
            & threshold_in=solver_parameters%threshold)
       CALL TransposeMatrix(StackV, StackVT)
       IF (StackVT%is_complex) THEN
          CALL ConjugateMatrix(StackVT)
       END IF
       CALL MatrixMultiply(StackVT, TempMat, VAV, &
            & threshold_in=solver_parameters%threshold)
       CALL DestructMatrix(StackVT)
       CALL DestructMatrix(TempMat)

       !! Iterate Recursively
       CALL ExtractCorner(VAV, left_dim, right_dim, LeftMat, RightMat)
       CALL DestructMatrix(VAV)
       CALL EigenRecursive(LeftMat, LeftVectors, splits(:midpoint-1), &
            & solver_parameters)
       CALL DestructMatrix(LeftMat)
       CALL EigenRecursive(RightMat, RightVectors, &
            & splits(midpoint+1:)-left_dim, solver_parameters)
       CALL DestructMatrix(RightMat)
       CALL ConstructEmptyMatrix(TempMat, this)
       CALL StackMatrices(LeftVectors, RightVectors, left_dim, left_dim, &
            & TempMat)
       CALL DestructMatrix(LeftVectors)
       CALL DestructMatrix(RightVectors)

       !! Recombine
       CALL MatrixMultiply(StackV, TempMat, eigenvectors, &
            & threshold_in=solver_parameters%threshold)

       !! Cleanup
       CALL DestructMatrix(StackV)
       CALL DestructMatrix(TempMat)
    END IF
    CALL DestructMatrix(Identity)

  END SUBROUTINE EigenRecursive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given two matrices, this stacks one on the left and right. The offsets
  !> describe how to shift the values of the right mat.
  SUBROUTINE StackMatrices(LeftMat,RightMat,column_offset,row_offset,Combined)
    !> The matrix on the left.
    TYPE(Matrix_ps), INTENT(IN) :: LeftMat
    !> The matrix on the right.
    TYPE(Matrix_ps), INTENT(IN) :: RightMat
    !> How many columns over to shift
    INTEGER :: column_offset
    !> How many rows over to shift
    INTEGER :: row_offset
    !> The stacked matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: Combined
    !! Local Data
    TYPE(TripletList_r) :: left_triplets, right_triplets, combined_triplets
    INTEGER :: left_size, right_size, total_size
    INTEGER :: counter

    !! Basic Triplet Lists
    CALL GetMatrixTripletList(LeftMat, left_triplets)
    CALL GetMatrixTripletList(RightMat, right_triplets)
    left_size = left_triplets%CurrentSize
    right_size = right_triplets%CurrentSize
    total_size = left_size + right_size
    combined_triplets = TripletList_r(total_size)

    !! Adjust right triplets
    DO counter = 1, right_size
       right_triplets%data(counter)%index_column = &
            & right_triplets%data(counter)%index_column  + column_offset
       right_triplets%data(counter)%index_row = &
            & right_triplets%data(counter)%index_row  + row_offset
    END DO

    !! Combine
    combined_triplets%data(:left_size) = left_triplets%data(:left_size)
    combined_triplets%data(left_size+1:total_size) = &
         & right_triplets%data(:right_size)
    CALL FillMatrixFromTripletList(Combined, combined_triplets, .TRUE.)

    !! Cleanup
    CALL DestructTripletList(left_triplets)
    CALL DestructTripletList(right_triplets)
    CALL DestructTripletList(combined_triplets)

  END SUBROUTINE StackMatrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract the top left and the bottom right corners of a matrix.
  SUBROUTINE ExtractCorner(this, left_dim, right_dim, LeftMat, RightMat)
    !> The matrix to extract from.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The top left corner.
    TYPE(Matrix_ps), INTENT(INOUT) :: LeftMat
    !> The bottom right corner.
    TYPE(Matrix_ps), INTENT(INOUT) :: RightMat
    !> The dimension of the left corner.
    INTEGER :: left_dim
    !> The dimension of the right corner.
    INTEGER :: right_dim
    !! Local Data
    TYPE(TripletList_r) :: left_triplets, right_triplets, combined_triplets
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: counter

    !! Get Triplet Lists
    CALL GetMatrixTripletList(this, combined_triplets)
    left_triplets = TripletList_r()
    right_triplets = TripletList_r()

    !! Extract Triplets
    DO counter = 1, combined_triplets%CurrentSize
       CALL GetTripletAt(combined_triplets, counter, temp_triplet)
       IF (temp_triplet%index_row .LE. left_dim .AND. &
            & temp_triplet%index_column .LE. left_dim) THEN
          CALL AppendToTripletList(left_triplets, temp_triplet)
       ELSE IF (temp_triplet%index_row .GT. left_dim .AND. &
            & temp_triplet%index_column .GT. left_dim) THEN
          temp_triplet%index_row = temp_triplet%index_row - left_dim
          temp_triplet%index_column = temp_triplet%index_column - left_dim
          CALL AppendToTripletList(right_triplets, temp_triplet)
       END IF
    END DO

    !! Fill
    CALL ConstructEmptyMatrix(LeftMat, left_dim, this%process_grid, &
         & this%is_complex)
    CALL ConstructEmptyMatrix(RightMat, right_dim, this%process_grid, &
         & this%is_complex)
    CALL FillMatrixFromTripletList(LeftMat, left_triplets, .TRUE.)
    CALL FillMatrixFromTripletList(RightMat, right_triplets, .TRUE.)

    !! Cleanup
    CALL DestructTripletList(left_triplets)
    CALL DestructTripletList(combined_triplets)
    CALL DestructTripletList(right_triplets)

  END SUBROUTINE ExtractCorner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> When we want to only compute the first n eigenvalues of a matrix, this
  !> routine will project out the higher eigenvalues.
  SUBROUTINE ReduceDimension(this, dim, parameters, ReducedMat)
    !> The starting matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The number of eigenvalues ot keep.
    INTEGER, INTENT(IN) :: dim
    !> a dimxdim matrix with the same first n eigenvalues as the first.
    TYPE(Matrix_ps), INTENT(INOUT) :: ReducedMat
    !> The solver parameters.
    TYPE(SolverParameters_t), INTENT(IN) :: parameters
    !! Local Variables - matrices
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: PMat
    TYPE(Matrix_ps) :: PVec, PVecT
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: VAV
    !! Special parameters
    TYPE(SolverParameters_t) :: balanced_it
    TYPE(Permutation_t) :: permutation

    !! Compute Identity Matrix
    CALL ConstructEmptyMatrix(Identity, this)
    CALL FillMatrixIdentity(Identity)

    !! Setup special parameters for purification
    CALL ConstructRandomPermutation(permutation, &
         & this%logical_matrix_dimension)
    balanced_it = SolverParameters_t(&
         & converge_diff_in = parameters%converge_diff, &
         & threshold_in = parameters%threshold, &
         & max_iterations_in = parameters%max_iterations, &
         & BalancePermutation_in = permutation)

    !! Purify
    CALL TRS2(this, Identity, dim*2, PMat, solver_parameters_in=balanced_it)

    !! Compute Eigenvectors of the Density Matrix
    CALL PivotedCholeskyDecomposition(PMat, PVec, dim, &
         & solver_parameters_in=parameters)

    !! Rotate to the divided subspace
    CALL MatrixMultiply(this, PVec, TempMat, threshold_in=parameters%threshold)
    CALL TransposeMatrix(PVec, PVecT)
    IF (PVecT%is_complex) THEN
       CALL ConjugateMatrix(PVecT)
    END IF
    CALL MatrixMultiply(PVecT, TempMat, VAV, threshold_in=parameters%threshold)

    !! Extract
    CALL ExtractCorner(VAV, dim, dim, ReducedMat, TempMat)

    CALL DestructMatrix(Identity)
    CALL DestructMatrix(PMat)
    CALL DestructMatrix(PVec)
    CALL DestructMatrix(PVecT)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(VAV)
    CALL DestructPermutation(permutation)
  END SUBROUTINE ReduceDimension
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
  !> Figure out where to split the matrix
  SUBROUTINE FillSplits(this, splits)
    !> The matrix we are computing.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The splits are recorded in this vector.
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: splits
    !! Local Variables
    TYPE(Matrix_ps) :: vec, val
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: val_array
    TYPE(TripletList_r) :: diag_list
    INTEGER :: II
    INTEGER :: split_count
    REAL(NTREAL) :: split_min = 0.01
    REAL(NTREAL) :: range
    INTEGER :: ierr

    !! "Guess" The Eigenvalues - test
    CALL ReferenceEigenDecomposition(this, vec, val)
    CALL GetMatrixTripletList(val, diag_list)
    ALLOCATE(val_array(this%actual_matrix_dimension))
    val_array = 0
    DO II = 1, diag_list%CurrentSize
       val_array(diag_list%data(II)%index_row) = diag_list%data(II)%point_value
    END DO
    CALL MPI_Allreduce(MPI_IN_PLACE, val_array, this%actual_matrix_dimension, &
         & MPINTREAL, MPI_SUM, this%process_grid%within_slice_comm, ierr)
    range = val_array(SIZE(val_array))-val_array(1)

    !! Fill The Split Array
    split_count = 0
    DO II = 1, SIZE(val_array)-1
       IF (ABS(val_array(II+1) - val_array(II))/range .GT. split_min) THEN
          split_count = split_count + 1
       END IF
    END DO
    ALLOCATE(splits(split_count))
    split_count = 1
    DO II = 1, SIZE(val_array)-1
       IF (ABS(val_array(II+1) - val_array(II))/range .GT. split_min) THEN
          splits(split_count) = II
          split_count = split_count + 1
       END IF
    END DO

    !! Cleanup
    CALL DestructTripletList(diag_list)
    DEALLOCATE(val_array)

  END SUBROUTINE FillSplits
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
    TYPE(TripletList_r) :: triplet_list, sorted_triplet_list, triplet_w
    TYPE(TripletList_r), DIMENSION(:), ALLOCATABLE :: send_list
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
    TYPE(TripletList_c) :: triplet_list, sorted_triplet_list, triplet_w
    TYPE(TripletList_c), DIMENSION(:), ALLOCATABLE :: send_list
    TYPE(Matrix_lsc) :: sparse
    TYPE(Matrix_lsc) :: local_a, local_v
    TYPE(Matrix_ldc) :: dense_a, dense_v, dense_w

    INCLUDE "solver_includes/EigenSerial.f90"
  END SUBROUTINE EigenSerial_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule
