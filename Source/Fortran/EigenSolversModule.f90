!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the eigenvalues or singular values of a matrix.
MODULE EigenSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE DenseMatrixModule, ONLY : DenseMatrix_t, &
       & ConstructSparseFromDense, ConstructDenseFromSparse, &
       & DenseEigenDecomposition, DestructDenseMatrix
  USE DistributedSparseMatrixAlgebraModule, ONLY : DistributedGemm, &
       & IncrementDistributedSparseMatrix
  USE DistributedSparseMatrixModule, ONLY : DistributedSparseMatrix_t, &
       & CommSplitDistributedSparseMatrix, &
       & ConstructEmptyDistributedSparseMatrix, CopyDistributedSparseMatrix, &
       & DestructDistributedSparseMatrix, FillFromTripletList, &
       & FillDistributedIdentity, GetTripletList, PrintMatrixInformation, &
       & TransposeDistributedSparseMatrix, GetSize
  USE DensityMatrixSolversModule, ONLY : TRS2
  USE FixedSolversModule, ONLY : FixedSolverParameters_t
  USE IterativeSolversModule, ONLY : IterativeSolverParameters_t, &
       PrintIterativeSolverParameters
#if EIGENEXA
  USE EigenExaModule, ONLY : EigenExa_s
#endif
  USE LinearSolversModule, ONLY : PivotedCholeskyDecomposition
  USE PermutationModule, ONLY : ConstructRandomPermutation, Permutation_t
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteElement, WriteListElement
  USE ParameterConverterModule, ONLY : ConvertIterativeToFixed
  USE SignSolversModule, ONLY : PolarDecomposition
  USE SparseMatrixModule, ONLY : SparseMatrix_t, ConstructFromTripletList, &
       & DestructSparseMatrix, MatrixToTripletList
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY : TripletList_t, AppendToTripletList, &
       & ConstructTripletList, DestructTripletList, GetTripletAt, &
       & RedistributeTripletLists, SortTripletList
  USE TripletModule, ONLY : Triplet_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, PARAMETER :: BASESIZE = 2048
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReferenceEigenDecomposition
  PUBLIC :: SingularValueDecomposition
  PUBLIC :: SplittingEigenDecomposition
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine uses a dense eigensolver for reference purposes.
  !! @param[in] this the matrix to decompose.
  !! @param[inout] eigenvectors the eigenvectors of a matrix.
  !! @param[in] solver_parameters_in parameters for computing (optional).
  SUBROUTINE ReferenceEigenDecomposition(this, eigenvectors, &
       & solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! For Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: iter_params
    TYPE(FixedSolverParameters_t) :: fixed_params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       iter_params = solver_parameters_in
    ELSE
       iter_params = IterativeSolverParameters_t()
    END IF
    CALL ConvertIterativeToFixed(iter_params, fixed_params)

#if EIGENEXA
    CALL EigenExa_s(this, eigenvectors, fixed_params)
#else
    CALL SerialBase(this, eigenvectors, fixed_params)
#endif

  END SUBROUTINE ReferenceEigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvalues and eigenvectors of a matrix.
  !! @param[in] this the matrix to decompose.
  !! @param[out] eigenvectors a matrix containing the eigenvectors of a matrix.
  !! @param[out] eigenvalues a diagonal matrix containing the eigenvalues.
  !! @param[in] num_values_in the number of eigenvalues to compute (optional).
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE SplittingEigenDecomposition(this, eigenvectors, eigenvalues, &
       & num_values_in, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvalues
    INTEGER, INTENT(IN), OPTIONAL :: num_values_in
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    TYPE(FixedSolverParameters_t) :: fixed_param
    TYPE(IterativeSolverParameters_t) :: it_param
    INTEGER :: num_values
    !! Local
    TYPE(DistributedSparseMatrix_t) :: eigenvectorsT, TempMat
    TYPE(DistributedSparseMatrix_t) :: ReducedMat

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
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
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Setup the solver parameters
    it_param = IterativeSolverParameters_t(&
         & converge_diff_in = solver_parameters%converge_diff, &
         & threshold_in = solver_parameters%threshold, &
         & max_iterations_in = solver_parameters%max_iterations)
    CALL ConvertIterativeToFixed(it_param, fixed_param)

    !! Rotate so that we are only computing the target number of values
    CALL StartTimer("Initial Purify")
    IF (num_values .LT. this%actual_matrix_dimension) THEN
       CALL ReduceDimension(this, num_values, solver_parameters, &
            & fixed_param, ReducedMat)
    ELSE
       CALL CopyDistributedSparseMatrix(this, ReducedMat)
    END IF
    CALL StopTimer("Initial Purify")

    !! Actual Solve
    CALL EigenRecursive(ReducedMat, eigenvectors, solver_parameters, it_param, &
         & fixed_param)
    CALL TransposeDistributedSparseMatrix(eigenvectors, eigenvectorsT)
    CALL DistributedGemm(eigenvectorsT, ReducedMat, TempMat, &
         & threshold_in=solver_parameters%threshold)
    CALL DistributedGemm(TempMat, eigenvectors, eigenvalues, &
         & threshold_in=solver_parameters%threshold)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructDistributedSparseMatrix(eigenvectorsT)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(ReducedMat)
  END SUBROUTINE SplittingEigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the singular values and singular vectors of a matrix.
  !! @param[in] this the matrix to decompose.
  !! @param[out] left_vectors a matrix containing the left singular vectors.
  !! @param[out] right_vectors a matrix containing the right singular vectors.
  !! @param[out] singularvalues a diagonal matrix containing the singularvalues.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE SingularValueDecomposition(this, left_vectors, &
       & right_vectors, singularvalues, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: left_vectors
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: right_vectors
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: singularvalues
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(DistributedSparseMatrix_t) :: UMat, HMat

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Singular Value Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Recursive")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! First compute the polar decomposition of the matrix.
    CALL PolarDecomposition(this, UMat, HMat, solver_parameters)

    !! Compute the eigen decomposition of the hermitian matrix
    CALL SplittingEigenDecomposition(HMat, right_vectors, singularvalues, &
         & solver_parameters_in=solver_parameters)

    !! Compute the left singular vectors
    CALL DistributedGemm(UMat, right_vectors, left_vectors, &
         & threshold_in=solver_parameters%threshold)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructDistributedSparseMatrix(UMat)
    CALL DestructDistributedSparseMatrix(HMat)
  END SUBROUTINE SingularValueDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The recursive workhorse routine for the eigendecompositon.
  !! @param[in] this the matrix to decompose.
  !! @param[out] eigenvectors a matrix containing the eigenvectors of a matrix.
  !! @param[in] it_parameters parameters for the iterative solvers.
  !! @param[in] it_parameters parameters for the fixed solvers.
  RECURSIVE SUBROUTINE EigenRecursive(this, eigenvectors, solver_parameters, &
       & it_param, fixed_param)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(IterativeSolverParameters_t), INTENT(IN) :: solver_parameters
    TYPE(IterativeSolverParameters_t), INTENT(IN) :: it_param
    TYPE(FixedSolverParameters_t), INTENT(IN) :: fixed_param
    !! Local Variables - matrices
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: PMat, PHoleMat
    TYPE(DistributedSparseMatrix_t) :: PVec, PHoleVec
    TYPE(DistributedSparseMatrix_t) :: StackV, StackVT
    TYPE(DistributedSparseMatrix_t) :: VAV
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedSparseMatrix_t) :: LeftMat, RightMat
    TYPE(DistributedSparseMatrix_t) :: LeftVectors, RightVectors
    TYPE(DistributedSparseMatrix_t) :: SubMat, SubVec
    !! Local Variables - For Splitting
    INTEGER :: mat_dim, left_dim, right_dim
    INTEGER :: color
    LOGICAL :: split_slice
    REAL(NTREAL) :: sparsity
    !! Special parameters
    TYPE(IterativeSolverParameters_t) :: balanced_it
    TYPE(Permutation_t) :: permutation

    mat_dim = this%actual_matrix_dimension
    CALL ConstructEmptyDistributedSparseMatrix(Identity, mat_dim, &
         & process_grid_in=this%process_grid)
    CALL FillDistributedIdentity(Identity)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteListElement(key="Iteration_Size", int_value_in=mat_dim)
       CALL EnterSubLog
       CALL PrintMatrixInformation(this)
       CALL ExitSubLog
    END IF
    sparsity = REAL(GetSize(this),KIND=NTREAL) / &
         & (REAL(this%actual_matrix_dimension,KIND=NTREAL)**2)

    !! Base Case
    IF (mat_dim .LE. BASESIZE .OR. sparsity .GT. 0.30) THEN
       CALL StartTimer("Base Case")
#if EIGENEXA
       CALL EigenExa_s(this, eigenvectors, fixed_param)
#else
       CALL SerialBase(this, eigenvectors, fixed_param)
#endif
       CALL StopTimer("Base Case")
    ELSE
       !! Setup
       left_dim = mat_dim/2
       right_dim = mat_dim - left_dim

       !! Setup special parameters for purification
       CALL ConstructRandomPermutation(permutation, &
            & this%logical_matrix_dimension)
       balanced_it = IterativeSolverParameters_t(&
            & converge_diff_in = it_param%converge_diff, &
            & threshold_in = it_param%threshold, &
            & max_iterations_in = it_param%max_iterations, &
            & BalancePermutation_in=permutation)

       !! Purify
       CALL StartTimer("Purify")
       CALL TRS2(this,Identity,left_dim*2,PMat,solver_parameters_in=balanced_it)
       CALL CopyDistributedSparseMatrix(Identity, PHoleMat)
       CALL IncrementDistributedSparseMatrix(PMat, PHoleMat, &
            & alpha_in=REAL(-1.0,NTREAL), threshold_in=it_param%threshold)
       CALL StopTimer("Purify")

       !! Compute Eigenvectors of the Density Matrix
       CALL StartTimer("Cholesky")
       CALL PivotedCholeskyDecomposition(PMat, PVec, left_dim, &
            & solver_parameters_in=fixed_param)
       CALL PivotedCholeskyDecomposition(PHoleMat, PHoleVec, right_dim, &
            & solver_parameters_in=fixed_param)
       CALL StopTimer("Cholesky")
       CALL ConstructEmptyDistributedSparseMatrix(StackV, &
            & this%actual_matrix_dimension, process_grid_in=this%process_grid)
       CALL StartTimer("Stack")
       CALL StackMatrices(PVec, PHoleVec, left_dim, 0, StackV)
       CALL StopTimer("Stack")

       !! Rotate to the divided subspace
       CALL StartTimer("Rotate")
       CALL DistributedGemm(this, StackV, TempMat, &
            & threshold_in=it_param%threshold)
       CALL TransposeDistributedSparseMatrix(StackV, StackVT)
       CALL DistributedGemm(StackVT, TempMat, VAV, &
            & threshold_in=it_param%threshold)
       CALL StopTimer("Rotate")

       !! Iterate Recursively
       IF (this%process_grid%total_processors .GT. 1) THEN
          CALL ExtractCornerS(VAV, left_dim, right_dim, SubMat, color, &
               & split_slice)
          CALL EigenRecursive(SubMat, SubVec, solver_parameters, it_param,&
               & fixed_param)
          CALL ConstructEmptyDistributedSparseMatrix(TempMat, &
               & this%actual_matrix_dimension, this%process_grid)
          CALL StackMatricesS(SubVec, left_dim, left_dim, color, TempMat)
       ELSE
          CALL ExtractCorner(VAV, left_dim, right_dim, LeftMat, RightMat)
          CALL EigenRecursive(LeftMat,LeftVectors,solver_parameters, it_param,&
               & fixed_param)
          CALL EigenRecursive(RightMat,RightVectors,solver_parameters,it_param,&
               & fixed_param)
          CALL ConstructEmptyDistributedSparseMatrix(TempMat, &
               & this%actual_matrix_dimension, this%process_grid)
          CALL StackMatrices(LeftVectors, RightVectors, left_dim, left_dim, &
               & TempMat)
       END IF

       !! Recombine
       CALL StartTimer("Recombine")
       CALL DistributedGemm(StackV, TempMat, eigenvectors, &
            & threshold_in=it_param%threshold)
       CALL StopTimer("Recombine")

       !! Cleanup
       CALL DestructDistributedSparseMatrix(SubMat)
       CALL DestructDistributedSparseMatrix(SubVec)
       CALL DestructDistributedSparseMatrix(PMat)
       CALL DestructDistributedSparseMatrix(PHoleMat)
       CALL DestructDistributedSparseMatrix(PVec)
       CALL DestructDistributedSparseMatrix(PHoleVec)
       CALL DestructDistributedSparseMatrix(StackV)
       CALL DestructDistributedSparseMatrix(StackVT)
       CALL DestructDistributedSparseMatrix(VAV)
       CALL DestructDistributedSparseMatrix(TempMat)
       CALL DestructDistributedSparseMatrix(LeftMat)
       CALL DestructDistributedSparseMatrix(RightMat)
       CALL DestructDistributedSparseMatrix(LeftVectors)
       CALL DestructDistributedSparseMatrix(RightVectors)
    END IF
    CALL DestructDistributedSparseMatrix(Identity)

  END SUBROUTINE EigenRecursive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given two matrices, this stacks one on the left and right. The offsets
  !! describe how to shift the values of the right mat.
  !! @param[in] LeftMat the matrix on the left.
  !! @param[in] RightMat the matrix on the right.
  !! @param[in] column_offset how many columns over to shift
  !! @param[in] row_offset how many rows over to shift
  !! @param[out] Combined the stacked matrix.
  SUBROUTINE StackMatrices(LeftMat,RightMat,column_offset,row_offset,Combined)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: LeftMat, RightMat
    INTEGER :: column_offset, row_offset
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: Combined
    !! Local Data
    TYPE(TripletList_t) :: left_triplets, right_triplets, combined_triplets
    INTEGER :: left_size, right_size, total_size
    INTEGER :: counter

    !! Basic Triplet Lists
    CALL GetTripletList(LeftMat, left_triplets)
    CALL GetTripletList(RightMat, right_triplets)
    left_size = left_triplets%CurrentSize
    right_size = right_triplets%CurrentSize
    total_size = left_size + right_size
    CALL ConstructTripletList(combined_triplets, total_size)

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
    CALL FillFromTripletList(Combined, combined_triplets, .TRUE.)

    !! Cleanup
    CALL DestructTripletList(left_triplets)
    CALL DestructTripletList(right_triplets)
    CALL DestructTripletList(combined_triplets)

  END SUBROUTINE StackMatrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE StackMatricesS(SubMat, column_offset, row_offset, color, FullMat)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: SubMat
    INTEGER, INTENT(IN) :: column_offset, row_offset
    INTEGER, INTENT(IN) :: color
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: FullMat
    !! Local Variables
    TYPE(TripletList_t) :: sub_triplets
    INTEGER :: counter

    !! Basic Triplet Lists
    CALL GetTripletList(SubMat, sub_triplets)

    !! Adjust right triplets
    IF (color .EQ. 1) THEN
       DO counter = 1, sub_triplets%CurrentSize
          sub_triplets%data(counter)%index_column = &
               & sub_triplets%data(counter)%index_column  + column_offset
          sub_triplets%data(counter)%index_row = &
               & sub_triplets%data(counter)%index_row  + row_offset
       END DO
    END IF
    IF (SubMat%process_grid%my_slice .GT. 0) THEN
       CALL ConstructTripletList(sub_triplets)
    END IF

    !! Combine
    CALL FillFromTripletList(FullMat, sub_triplets, .FALSE.)

    !! Cleanup
    CALL DestructTripletList(sub_triplets)

  END SUBROUTINE StackMatricesS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract the top left and the bottom right corners of a matrix.
  !! @param[in] this the matrix to extract from.
  !! @param[in] left_dim the dimension of the left corner.
  !! @param[in] right_dim the dimension of the right corner.
  !! @param[out] LeftMat the top left corner matrix.
  !! @param[out] RightMat the bottom right corner matirx.
  SUBROUTINE ExtractCorner(this, left_dim, right_dim, LeftMat, RightMat)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: LeftMat, RightMat
    INTEGER :: left_dim, right_dim
    !! Local Data
    TYPE(TripletList_t) :: left_triplets, right_triplets, combined_triplets
    TYPE(Triplet_t) :: temp_triplet
    INTEGER :: counter

    !! Get Triplet Lists
    CALL GetTripletList(this, combined_triplets)
    CALL ConstructTripletList(left_triplets)
    CALL ConstructTripletList(right_triplets)

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
    CALL ConstructEmptyDistributedSparseMatrix(LeftMat, left_dim, &
         & this%process_grid)
    CALL ConstructEmptyDistributedSparseMatrix(RightMat, right_dim, &
         & this%process_grid)
    CALL FillFromTripletList(LeftMat, left_triplets, .TRUE.)
    CALL FillFromTripletList(RightMat, right_triplets, .TRUE.)

    !! Cleanup
    CALL DestructTripletList(left_triplets)
    CALL DestructTripletList(combined_triplets)
    CALL DestructTripletList(right_triplets)

  END SUBROUTINE ExtractCorner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ExtractCornerS(this, left_dim, right_dim, SubMat, color, &
       & split_slice)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: SubMat
    LOGICAL, INTENT(OUT) :: split_slice
    INTEGER, INTENT(IN) :: left_dim, right_dim
    INTEGER, INTENT(OUT) :: color
    !! Local Data
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(TripletList_t) :: full_triplets, extracted_triplets
    TYPE(Triplet_t) :: temp_triplet
    INTEGER :: counter

    !! First Duplicate Across Process Grids
    CALL CommSplitDistributedSparseMatrix(this, TempMat, color, split_slice)

    !! Extract The Corner
    CALL GetTripletList(TempMat, full_triplets)

    !! Extract Triplets
    CALL ConstructTripletList(extracted_triplets)
    DO counter = 1, full_triplets%CurrentSize
       CALL GetTripletAt(full_triplets, counter, temp_triplet)
       IF (temp_triplet%index_row .LE. left_dim .AND. &
            & temp_triplet%index_column .LE. left_dim .AND. color .EQ. 0) THEN
          CALL AppendToTripletList(extracted_triplets, temp_triplet)
       ELSE IF (temp_triplet%index_row .GT. left_dim .AND. &
            & temp_triplet%index_column .GT. left_dim .AND. color .EQ. 1) THEN
          temp_triplet%index_row = temp_triplet%index_row - left_dim
          temp_triplet%index_column = temp_triplet%index_column - left_dim
          CALL AppendToTripletList(extracted_triplets, temp_triplet)
       END IF
    END DO

    !! Fill
    IF (color .EQ. 0) THEN
       CALL ConstructEmptyDistributedSparseMatrix(SubMat, left_dim, &
            & process_grid_in=TempMat%process_grid)
    ELSE
       CALL ConstructEmptyDistributedSparseMatrix(SubMat, right_dim, &
            & process_grid_in=TempMat%process_grid)
    END IF
    CALL FillFromTripletList(SubMat,extracted_triplets,preduplicated_in=.TRUE.)

    !! Cleanup
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructTripletList(full_triplets)
    CALL DestructTripletList(extracted_triplets)

  END SUBROUTINE ExtractCornerS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> When we want to only compute the first n eigenvalues of a matrix, this
  !! routine will project out the higher eigenvalues.
  !! @param[in] this the starting matrix.
  !! @param[in] dim the number of eigenvalues ot keep.
  !! @param[in] it_param the iterative solver parameters.
  !! @param[in] fixed_param the fixed solver parameters.
  !! @param[out] ReducedMat a dimxdim matrix with the same first n eigenvalues
  !! as the first.
  SUBROUTINE ReduceDimension(this, dim, it_param, fixed_param, ReducedMat)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: dim
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: ReducedMat
    TYPE(IterativeSolverParameters_t), INTENT(IN) :: it_param
    TYPE(FixedSolverParameters_t), INTENT(IN) :: fixed_param
    !! Local Variables - matrices
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: PMat
    TYPE(DistributedSparseMatrix_t) :: PVec, PVecT
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedSparseMatrix_t) :: VAV
    !! Special parameters
    TYPE(IterativeSolverParameters_t) :: balanced_it
    TYPE(Permutation_t) :: permutation

    !! Compute Identity Matrix
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & this%actual_matrix_dimension, process_grid_in=this%process_grid)
    CALL FillDistributedIdentity(Identity)

    !! Setup special parameters for purification
    CALL ConstructRandomPermutation(permutation, &
         & this%logical_matrix_dimension)
    balanced_it = IterativeSolverParameters_t(&
         & converge_diff_in = it_param%converge_diff, &
         & threshold_in = it_param%threshold, &
         & max_iterations_in = it_param%max_iterations, &
         & BalancePermutation_in=permutation)

    !! Purify
    CALL TRS2(this, Identity, dim*2, PMat, solver_parameters_in=balanced_it)

    !! Compute Eigenvectors of the Density Matrix
    CALL PivotedCholeskyDecomposition(PMat, PVec, dim, &
         & solver_parameters_in=fixed_param)

    !! Rotate to the divided subspace
    CALL DistributedGemm(this, PVec, TempMat, &
         & threshold_in=it_param%threshold)
    CALL TransposeDistributedSparseMatrix(PVec, PVecT)
    CALL DistributedGemm(PVecT, TempMat, VAV, &
         & threshold_in=it_param%threshold)

    !! Extract
    CALL ExtractCorner(VAV, dim, dim, ReducedMat, TempMat)

    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(PMat)
    CALL DestructDistributedSparseMatrix(PVec)
    CALL DestructDistributedSparseMatrix(PVecT)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(VAV)
  END SUBROUTINE ReduceDimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve
  !! @param[in] this the matrix to compute
  !! @param[in] fixed_param the solve parameters.
  !! @param[out] eigenvecotrs the eigenvectors of the matrix
  SUBROUTINE SerialBase(this, eigenvectors, fixed_param)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(FixedSolverParameters_t), INTENT(IN) :: fixed_param
    !! Local Data
    TYPE(TripletList_t) :: triplet_list, sorted_triplet_list
    TYPE(TripletList_t), DIMENSION(:), ALLOCATABLE :: send_list
    TYPE(SparseMatrix_t) :: sparse
    INTEGER :: counter, list_size
    INTEGER :: mat_dim
    TYPE(SparseMatrix_t) :: local_a, local_v
    TYPE(DenseMatrix_t) :: dense_a, dense_v

    mat_dim = this%actual_matrix_dimension

    !! Gather on a single processor
    CALL GetTripletList(this, triplet_list)
    ALLOCATE(send_list(this%process_grid%slice_size))
    CALL ConstructTripletList(send_list(1), triplet_list%CurrentSize)
    DO counter = 2, this%process_grid%slice_size
       CALL ConstructTripletList(send_list(counter))
    END DO
    list_size = triplet_list%CurrentSize
    send_list(1)%data(:list_size) = triplet_list%data(:list_size)
    CALL DestructTripletList(triplet_list)
    CALL RedistributeTripletLists(send_list, &
         & this%process_grid%within_slice_comm, triplet_list)

    IF (this%process_grid%within_slice_rank .EQ. 0) THEN
       CALL SortTripletList(triplet_list, mat_dim, sorted_triplet_list, .TRUE.)
       CALL ConstructFromTripletList(local_a, sorted_triplet_list, mat_dim, &
            & mat_dim)
       CALL ConstructDenseFromSparse(local_a, dense_a)
       CALL StartTimer("Serial Inner")
       CALL DenseEigenDecomposition(dense_a, dense_v)
       CALL StopTimer("Serial Inner")
       CALL ConstructSparseFromDense(dense_v,local_v,fixed_param%threshold)
       CALL DestructDenseMatrix(dense_a)
       CALL DestructDenseMatrix(dense_v)
       CALL MatrixToTripletList(local_v,triplet_list)
    END IF
    CALL ConstructEmptyDistributedSparseMatrix(eigenvectors, mat_dim, &
         & process_grid_in=this%process_grid)
    CALL FillFromTripletList(eigenvectors, triplet_list, .TRUE.)

    !! Cleanup
    CALL DestructSparseMatrix(sparse)
    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
    DO counter = 1, this%process_grid%slice_size
       CALL DestructTripletList(send_list(counter))
    END DO
    DEALLOCATE(send_list)
  END SUBROUTINE SerialBase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule
