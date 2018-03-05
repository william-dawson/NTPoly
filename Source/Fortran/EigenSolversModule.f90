!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing the eigenvalues or singular values of a matrix.
MODULE EigenSolversModule
  USE DataTypesModule
  USE DistributedSparseMatrixAlgebraModule
  USE DistributedSparseMatrixModule
  USE DensityMatrixSolversModule
  USE FixedSolversModule
  USE IterativeSolversModule
  USE LinearSolversModule
  USE LoggingModule
  USE ParameterConverterModule
  USE PermutationModule
  USE TripletListModule
  USE TripletModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition
  ! PUBLIC :: SingularValueDecompostion
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvalues and eigenvectors of a matrix.
  SUBROUTINE EigenDecomposition(this, eigenvectors, eigenvalues, &
       & solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvalues
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    TYPE(FixedSolverParameters_t) :: f_solver_parameters
    TYPE(Permutation_t) :: default_perm
    !! Local
    TYPE(DistributedSparseMatrix_t) :: eigenvectorsT, TempMat

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Eigen Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="recursive")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Setup the solver parameters
    solver_parameters%be_verbose = .FALSE.
    CALL ConstructDefaultPermutation(default_perm,this%logical_matrix_dimension)
    CALL SetIterativeLoadBalance(solver_parameters, default_perm)
    CALL ConvertIterativeToFixed(solver_parameters, f_solver_parameters)

    !! Actual Solve
    CALL EigenRecursive(this, eigenvectors, solver_parameters, &
         & f_solver_parameters)
    CALL TransposeDistributedSparseMatrix(eigenvectors, eigenvectorsT)
    CALL DistributedGemm(eigenvectorsT, this, TempMat, &
         & threshold_in=solver_parameters%threshold)
    CALL DistributedGemm(TempMat, eigenvectors, eigenvalues, &
         & threshold_in=solver_parameters%threshold)

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructDistributedSparseMatrix(eigenvectorsT)
    CALL DestructDistributedSparseMatrix(TempMat)
  END SUBROUTINE EigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE EigenRecursive(this, eigenvectors, it_param, fixed_param)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
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
    !! Local Variables - For Splitting
    INTEGER :: mat_dim, left_dim, right_dim

    mat_dim = this%actual_matrix_dimension
    CALL ConstructEmptyDistributedSparseMatrix(Identity, mat_dim)
    CALL FillDistributedIdentity(Identity)

    !! Base Case
    IF (mat_dim .EQ. 1) THEN
       CALL CopyDistributedSparseMatrix(Identity, eigenvectors)
    ELSE
       !! Setup
       left_dim = mat_dim/2
       right_dim = mat_dim - left_dim

       !! Purify
       CALL TRS2(this,Identity,left_dim*2,PMat,solver_parameters_in=it_param)
       CALL CopyDistributedSparseMatrix(Identity, PHoleMat)
       CALL IncrementDistributedSparseMatrix(PMat, PHoleMat, &
            & alpha_in=REAL(-1.0,NTREAL), threshold_in=it_param%threshold)

       !! Compute Eigenvectors of the Density Matrix
       CALL PivotedCholeskyDecomposition(PMat, PVec, left_dim, &
            & solver_parameters_in=fixed_param)
       CALL PivotedCholeskyDecomposition(PHoleMat, PHoleVec, right_dim, &
            & solver_parameters_in=fixed_param)
       CALL ConstructEmptyDistributedSparseMatrix(StackV, &
            & this%actual_matrix_dimension)
       CALL StackMatrices(PVec, PHoleVec, left_dim, 0, StackV)

       !! Rotate to the divided subspace
       CALL DistributedGemm(this, StackV, TempMat, &
            & threshold_in=it_param%threshold)
       CALL TransposeDistributedSparseMatrix(StackV, StackVT)
       CALL DistributedGemm(StackVT, TempMat, VAV, &
            & threshold_in=it_param%threshold)

       !! Iterate Recursively
       CALL ExtractCorner(VAV, left_dim, right_dim, LeftMat, RightMat)
       CALL EigenRecursive(LeftMat,LeftVectors,it_param,fixed_param)
       CALL EigenRecursive(RightMat,RightVectors,it_param,fixed_param)

       !! Recombine
       CALL ConstructEmptyDistributedSparseMatrix(TempMat, &
            & this%actual_matrix_dimension)
       CALL StackMatrices(LeftVectors, RightVectors, left_dim, left_dim, &
            & TempMat)
       CALL DistributedGemm(StackV, TempMat, eigenvectors, &
            & threshold_in=it_param%threshold)

       !! Cleanup
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
    CALL ConstructEmptyDistributedSparseMatrix(LeftMat, left_dim)
    CALL ConstructEmptyDistributedSparseMatrix(RightMat, right_dim)
    CALL FillFromTripletList(LeftMat, left_triplets, .TRUE.)
    CALL FillFromTripletList(RightMat, right_triplets, .TRUE.)

    !! Cleanup
    CALL DestructTripletList(left_triplets)
    CALL DestructTripletList(right_triplets)

  END SUBROUTINE ExtractCorner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenSolversModule
