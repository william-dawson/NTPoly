!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the eigenvalues or singular values of a matrix.
MODULE EigenSolversModule
  USE DataTypesModule
  USE DensityMatrixSolversModule
  USE DMatrixModule
#if EIGENEXA
  USE EigenExaModule, ONLY : EigenExa_s
#endif
  USE LinearSolversModule
  USE PermutationModule
  USE LoggingModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE SMatrixModule
  USE SMatrixAlgebraModule
  USE SignSolversModule
  USE TimerModule
  USE TripletListModule
  USE TripletModule
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
  SUBROUTINE ReferenceEigenDecomposition(this, eigenvectors, eigenvalues_in, &
       & solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_in
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
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
  !! @param[in] this the matrix to decompose.
  !! @param[out] eigenvectors a matrix containing the eigenvectors of a matrix.
  !! @param[out] eigenvalues a diagonal matrix containing the eigenvalues.
  !! @param[in] num_values_in the number of eigenvalues to compute (optional).
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE SplittingEigenDecomposition(this, eigenvectors, eigenvalues, &
       & num_values_in, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    INTEGER, INTENT(IN), OPTIONAL :: num_values_in
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    INTEGER :: num_values
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
       CALL PrintParameters(solver_parameters)
    END IF

    !! Rotate so that we are only computing the target number of values
    CALL StartTimer("Initial Purify")
    IF (num_values .LT. this%actual_matrix_dimension) THEN
       CALL ReduceDimension(this, num_values, solver_parameters, ReducedMat)
    ELSE
       CALL CopyMatrix(this, ReducedMat)
    END IF
    CALL StopTimer("Initial Purify")

    !! Actual Solve
    CALL EigenRecursive(ReducedMat, eigenvectors, solver_parameters)

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

    CALL DestructMatrix(eigenvectorsT)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(eigenvalues_r)
    CALL DestructMatrix(ReducedMat)
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
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: left_vectors
    TYPE(Matrix_ps), INTENT(INOUT) :: right_vectors
    TYPE(Matrix_ps), INTENT(INOUT) :: singularvalues
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
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
  !! @param[in] this the matrix to decompose.
  !! @param[out] eigenvectors a matrix containing the eigenvectors of a matrix.
  !! @param[in] solver_parameters for the solvers.
  RECURSIVE SUBROUTINE EigenRecursive(this, eigenvectors, solver_parameters)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
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
    TYPE(Matrix_ps) :: SubMat, SubVec
    !! Local Variables - For Splitting
    INTEGER :: mat_dim, left_dim, right_dim
    INTEGER :: color
    LOGICAL :: split_slice
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
    IF (mat_dim .LE. BASESIZE .OR. sparsity .GT. 0.30) THEN
       CALL StartTimer("Base Case")
       CALL ReferenceEigenDecomposition(this, eigenvectors, &
            & solver_parameters_in=solver_parameters)
       CALL StopTimer("Base Case")
    ELSE
       !! Setup
       left_dim = mat_dim/2
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
       CALL StartTimer("Purify")
       CALL TRS2(this,Identity,left_dim*2,PMat,solver_parameters_in=balanced_it)
       CALL CopyMatrix(Identity, PHoleMat)
       CALL IncrementMatrix(PMat, PHoleMat, alpha_in=REAL(-1.0,NTREAL), &
            & threshold_in=solver_parameters%threshold)
       CALL StopTimer("Purify")

       !! Compute Eigenvectors of the Density Matrix
       CALL StartTimer("Cholesky")
       CALL PivotedCholeskyDecomposition(PMat, PVec, left_dim, &
            & solver_parameters_in=solver_parameters)
       CALL PivotedCholeskyDecomposition(PHoleMat, PHoleVec, right_dim, &
            & solver_parameters_in=solver_parameters)
       CALL StopTimer("Cholesky")
       CALL ConstructEmptyMatrix(StackV, this)
       CALL StartTimer("Stack")
       CALL StackMatrices(PVec, PHoleVec, left_dim, 0, StackV)
       CALL StopTimer("Stack")

       !! Rotate to the divided subspace
       CALL StartTimer("Rotate")
       CALL MatrixMultiply(this, StackV, TempMat, &
            & threshold_in=solver_parameters%threshold)
       CALL TransposeMatrix(StackV, StackVT)
       IF (StackVT%is_complex) THEN
          CALL ConjugateMatrix(StackVT)
       END IF
       CALL MatrixMultiply(StackVT, TempMat, VAV, &
            & threshold_in=solver_parameters%threshold)
       CALL StopTimer("Rotate")

       !! Iterate Recursively
       IF (this%process_grid%total_processors .GT. 1 .AND. .FALSE.) THEN
          CALL ExtractCornerS(VAV, left_dim, right_dim, SubMat, color, &
               & split_slice)
          CALL EigenRecursive(SubMat, SubVec, solver_parameters)
          CALL ConstructEmptyMatrix(TempMat, this)
          CALL StackMatricesS(SubVec, left_dim, left_dim, color, TempMat)
       ELSE
          CALL ExtractCorner(VAV, left_dim, right_dim, LeftMat, RightMat)
          CALL EigenRecursive(LeftMat,LeftVectors,solver_parameters)
          CALL EigenRecursive(RightMat,RightVectors,solver_parameters)
          CALL ConstructEmptyMatrix(TempMat, this)
          CALL StackMatrices(LeftVectors, RightVectors, left_dim, left_dim, &
               & TempMat)
       END IF

       !! Recombine
       CALL StartTimer("Recombine")
       CALL MatrixMultiply(StackV, TempMat, eigenvectors, &
            & threshold_in=solver_parameters%threshold)
       CALL StopTimer("Recombine")

       !! Cleanup
       CALL DestructMatrix(SubMat)
       CALL DestructMatrix(SubVec)
       CALL DestructMatrix(PMat)
       CALL DestructMatrix(PHoleMat)
       CALL DestructMatrix(PVec)
       CALL DestructMatrix(PHoleVec)
       CALL DestructMatrix(StackV)
       CALL DestructMatrix(StackVT)
       CALL DestructMatrix(VAV)
       CALL DestructMatrix(TempMat)
       CALL DestructMatrix(LeftMat)
       CALL DestructMatrix(RightMat)
       CALL DestructMatrix(LeftVectors)
       CALL DestructMatrix(RightVectors)
    END IF
    CALL DestructMatrix(Identity)

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
    TYPE(Matrix_ps), INTENT(IN) :: LeftMat, RightMat
    INTEGER :: column_offset, row_offset
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
  SUBROUTINE StackMatricesS(SubMat, column_offset, row_offset, color, FullMat)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: SubMat
    INTEGER, INTENT(IN) :: column_offset, row_offset
    INTEGER, INTENT(IN) :: color
    TYPE(Matrix_ps), INTENT(INOUT) :: FullMat
    !! Local Variables
    TYPE(TripletList_r) :: sub_triplets
    INTEGER :: counter

    !! Basic Triplet Lists
    CALL GetMatrixTripletList(SubMat, sub_triplets)

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
       sub_triplets = TripletList_r()
    END IF

    !! Combine
    CALL FillMatrixFromTripletList(FullMat, sub_triplets, .FALSE.)

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
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: LeftMat, RightMat
    INTEGER :: left_dim, right_dim
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
  SUBROUTINE ExtractCornerS(this, left_dim, right_dim, SubMat, color, &
       & split_slice)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: SubMat
    LOGICAL, INTENT(OUT) :: split_slice
    INTEGER, INTENT(IN) :: left_dim, right_dim
    INTEGER, INTENT(OUT) :: color
    !! Local Data
    TYPE(Matrix_ps) :: TempMat
    TYPE(TripletList_r) :: full_triplets, extracted_triplets
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: counter

    !! First Duplicate Across Process Grids
    CALL CommSplitMatrix(this, TempMat, color, split_slice)

    !! Extract The Corner
    CALL GetMatrixTripletList(TempMat, full_triplets)

    !! Extract Triplets
    extracted_triplets = TripletList_r()
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
       CALL ConstructEmptyMatrix(SubMat, left_dim, TempMat%process_grid, &
            & TempMat%is_complex)
    ELSE
       CALL ConstructEmptyMatrix(SubMat, right_dim, TempMat%process_grid, &
            & TempMat%is_complex)
    END IF
    CALL FillMatrixFromTripletList(SubMat,extracted_triplets,preduplicated_in=.TRUE.)

    !! Cleanup
    CALL DestructMatrix(TempMat)
    CALL DestructTripletList(full_triplets)
    CALL DestructTripletList(extracted_triplets)

  END SUBROUTINE ExtractCornerS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> When we want to only compute the first n eigenvalues of a matrix, this
  !! routine will project out the higher eigenvalues.
  !! @param[in] this the starting matrix.
  !! @param[in] dim the number of eigenvalues ot keep.
  !! @param[in] parameters the solver parameters.
  !! @param[out] ReducedMat a dimxdim matrix with the same first n eigenvalues
  !! as the first.
  SUBROUTINE ReduceDimension(this, dim, parameters, ReducedMat)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: dim
    TYPE(Matrix_ps), INTENT(INOUT) :: ReducedMat
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
  END SUBROUTINE ReduceDimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve.
  !! @param[in] this the matrix to compute.
  !! @param[out] eigenvectors the eigenvectors of the matrix.
  !! @param[out] eigenvalues the eigenvalues of the matrix.
  !! @param[in] solver_parameters_in the solve parameters.
  SUBROUTINE EigenSerial(this, eigenvectors, eigenvalues_out, &
       & solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
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
  !> The base case: use lapack to solve.
  !! @param[in] this the matrix to compute.
  !! @param[out] eigenvectors the eigenvectors of the matrix.
  !! @param[out] eigenvalues the eigenvalues of the matrix.
  !! @param[in] solver_parameters_in the solve parameters.
  SUBROUTINE EigenSerial_r(this, eigenvectors, fixed_params, eigenvalues_out)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    TYPE(SolverParameters_t), INTENT(IN) :: fixed_params
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
    !! Local Data
    TYPE(TripletList_r) :: triplet_list, sorted_triplet_list, triplet_w
    TYPE(TripletList_r), DIMENSION(:), ALLOCATABLE :: send_list
    TYPE(Matrix_lsr) :: sparse
    TYPE(Matrix_lsr) :: local_a, local_v
    TYPE(Matrix_ldr) :: dense_a, dense_v, dense_w

    INCLUDE "SolverSupport/includes/EigenSerial.F90"
  END SUBROUTINE EigenSerial_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The base case: use lapack to solve.
  !! @param[in] this the matrix to compute.
  !! @param[out] eigenvectors the eigenvectors of the matrix.
  !! @param[out] eigenvalues the eigenvalues of the matrix.
  !! @param[in] solver_parameters_in the solve parameters.
  SUBROUTINE EigenSerial_c(this, eigenvectors, fixed_params, eigenvalues_out)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvectors
    TYPE(SolverParameters_t), INTENT(IN) :: fixed_params
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvalues_out
    !! Local Data
    TYPE(TripletList_c) :: triplet_list, sorted_triplet_list, triplet_w
    TYPE(TripletList_c), DIMENSION(:), ALLOCATABLE :: send_list
    TYPE(Matrix_lsc) :: sparse
    TYPE(Matrix_lsc) :: local_a, local_v
    TYPE(Matrix_ldc) :: dense_a, dense_v, dense_w

    INCLUDE "SolverSupport/includes/EigenSerial.F90"
  END SUBROUTINE EigenSerial_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule
