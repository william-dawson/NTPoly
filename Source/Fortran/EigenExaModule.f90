!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for calling eigenexa
MODULE EigenExaModule
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixModule, ONLY : DistributedSparseMatrix_t, &
       & GetTripletList, ConstructEmptyDistributedSparseMatrix, &
       & FillFromTripletList, PrintMatrixInformation
  USE FixedSolversModule, ONLY : FixedSolverParameters_t, &
       & PrintFixedSolverParameters
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader, WriteCitation
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletModule, ONLY : Triplet_t, SetTriplet
  USE TripletListModule, ONLY : TripletList_t, ConstructTripletList, &
       & GetTripletAt, DestructTripletList, AppendToTripletList, &
       & RedistributeTripletLists
  USE eigen_libs
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenExa_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, PRIVATE :: ExaHelper_t
     !> The number of processors involved.
     INTEGER :: num_procs
     !> The number of rows for the eigenexa comm.
     INTEGER :: proc_rows
     !> The number of columns for the eigenexa comm.
     INTEGER :: proc_cols
     !> Which process this is.
     INTEGER :: procid
     !> Which row is this process in.
     INTEGER :: rowid
     !> Which column is this process in.
     INTEGER :: colid
     !> The number of rows for the local matrix.
     INTEGER :: local_rows
     !> The number of columns for the local matrix.
     INTEGER :: local_cols
     !> The dimension fo the matrix.
     INTEGER :: mat_dim
     !> The communicator for this calculation.
     INTEGER :: comm
     !> Householder transform block size
     INTEGER :: MB
     !> Householder backward transformation block size
     INTEGER :: M
     !> Mode of the solver
     CHARACTER :: MODE
     INTEGER :: offset
  END TYPE ExaHelper_t
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a matrix using EigenExa.
  !! @param[in] A the matrix to decompose.
  !! @param[out] V the eigenvectors computed.
  !! @param[in] solver_parameters_in the parameters for this solver (optional).
  SUBROUTINE EigenExa_s(A, V, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: A
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: V
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Helper
    TYPE(ExaHelper_t) :: exa
    !! Dense Matrices
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: AD
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: VD
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: WD

    !! Process Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF

    !! Write info about the solver
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Eigen Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="EigenExa")
       CALL WriteCitation("imamura2011development")
       CALL PrintFixedSolverParameters(solver_parameters)
    END IF

    !! Main Routines
    CALL InitializeEigenExa(A, AD, VD, WD, exa)
    
    CALL StartTimer("NTToEigen")
    CALL NTToEigen(A, AD, exa)
    CALL StopTimer("NTToEigen")

    CALL StartTimer("EigenExaCompute")
    CALL Compute(AD, VD, WD, exa)
    CALL StopTimer("EigenExaCompute")

    CALL StartTimer("EigenToNT")
    CALL EigenToNT(VD, V, solver_parameters, exa)
    CALL StopTimer("EigenToNT")

    CALL CleanUp(AD, VD, WD)

    IF (solver_parameters%be_verbose) THEN
       CALL PrintMatrixInformation(V)
       CALL ExitSubLog
    END IF

  END SUBROUTINE EigenExa_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Setup the eigen exa data structures.
  !! @param[in] A the matrix we're working on.
  !! @param[inout] AD the dense matrix will be allocated.
  !! @param[inout] VD the dense matrix of eigenvectors will be allocated.
  !! @param[inout] WD the dense array of eigenvalues will be allocated.
  !! @param[inout] exa stores info about the calculation.
  SUBROUTINE InitializeEigenExa(A, AD, VD, WD, exa)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: A
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AD
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: VD
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WD
    TYPE(ExaHelper_t), INTENT(INOUT) :: exa
    !! Local Variables
    INTEGER :: ICTXT, INFO
    INTEGER, DIMENSION(9) :: DESCA
    INTEGER :: ierr

    !! Setup the MPI Communicator
    CALL MPI_Comm_dup(A%process_grid%global_comm, exa%comm, ierr)

    !! Build EigenExa Process Grid
    CALL eigen_init(exa%comm)
    CALL eigen_get_procs(exa%num_procs, exa%proc_rows, exa%proc_cols )
    CALL eigen_get_id(exa%procid, exa%rowid, exa%colid)

    !! Allocate Dense Matrices
    exa%mat_dim = A%actual_matrix_dimension
    CALL eigen_get_matdims(exa%mat_dim, exa%local_rows, exa%local_cols)
    ALLOCATE(AD(exa%local_rows, exa%local_cols))
    AD = 0
    ALLOCATE(VD(exa%local_rows, exa%local_cols))
    VD = 0
    ALLOCATE(WD(exa%mat_dim))
    WD = 0

    !> Default blocking parameters
    exa%MB = 128
    exa%M = 48
    exa%MODE = 'A'

    !! Blacs gives us the blocking info.
    ICTXT = eigen_get_blacs_context()
    CALL DESCINIT(DESCA, exa%mat_dim, exa%mat_dim, 1, 1, 0, 0, ICTXT, &
         & exa%local_rows, INFO )
    exa%offset = DESCA(9)

  END SUBROUTINE InitializeEigenExa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the distributed sparse matrix to a dense matrix in block-cyclic
  !! distribution.
  !! @param[in] A the matrix to convert.
  !! @param[inout] AD the dense, block-cyclic version.
  !! @param[inout] info about the calculation.
  SUBROUTINE NTToEigen(A, AD, exa)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: A
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: AD
    TYPE(ExaHelper_t), INTENT(INOUT) :: exa
    !! Local Variables
    TYPE(TripletList_t) :: triplet_a
    TYPE(TripletList_t), DIMENSION(:), ALLOCATABLE :: send_list
    TYPE(TripletList_t) :: recv_list
    TYPE(Triplet_t) :: trip, shifted_trip
    INTEGER :: ilookup, jlookup, iowner, jowner, ijowner
    INTEGER :: lrow, lcol
    INTEGER :: II

    !! We will fill a triplet list for each other process
    ALLOCATE(send_list(exa%num_procs))
    DO II = 1, exa%num_procs
       CALL ConstructTripletList(send_list(II))
    END DO

    !! Now Get The Triplet List, and adjust
    CALL GetTripletList(A, triplet_a)
    DO II = 1, triplet_a%CurrentSize
       CALL GetTripletAt(triplet_a, II, trip)

       !! Determine where that triplet will reside
       iowner = eigen_owner_node(trip%index_row, exa%proc_rows, exa%rowid)
       jowner = eigen_owner_node(trip%index_column, exa%proc_cols, exa%colid)
       ijowner = (jowner-1)*exa%proc_rows + iowner

       !! New indices
       ilookup = eigen_translate_g2l(trip%index_row, exa%proc_rows, &
            & exa%rowid)
       jlookup = eigen_translate_g2l(trip%index_column, exa%proc_cols, &
            & exa%colid)
       CALL SetTriplet(shifted_trip, jlookup, ilookup, trip%point_value)

       CALL AppendToTripletList(send_list(ijowner), shifted_trip)
    END DO

    !! Redistribute The Triplets
    CALL RedistributeTripletLists(send_list, exa%comm, recv_list)

    !! Write To The Dense Array
    DO II = 1, recv_list%CurrentSize
       lrow = recv_list%data(II)%index_row
       lcol = recv_list%data(II)%index_column
       AD(lrow,lcol) = recv_list%data(II)%point_value
    END DO

    !! Cleanup
    DO II = 1, exa%num_procs
       CALL DestructTripletList(send_list(II))
    END DO
    DEALLOCATE(send_list)
    CALL DestructTripletList(recv_list)

  END SUBROUTINE NTToEigen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the dense eigenvector matrix stored block-cyclicly back to
  !! a distributed sparse matrix.
  !! @param[in] VD the dense eigenvector matrix.
  !! @param[inout] V the distributed sparse matrix.
  !! @param[in] solver_parameters for thresholding small values.
  !! @param[inout] info about the calculation.
  SUBROUTINE EigenToNT(VD, V, solver_parameters, exa)
    !! Parameters
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: VD
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: V
    TYPE(FixedSolverParameters_t) :: solver_parameters
    TYPE(ExaHelper_t) :: exa
    !! Local Variables
    TYPE(TripletList_t) :: triplet_w, triplet_v
    TYPE(Triplet_t) :: trip
    INTEGER :: row_start, row_end, col_start, col_end
    INTEGER :: II, JJ, ilookup, jlookup
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: VD1
    INTEGER :: ind

    !! The Matrices We'll Build
    CALL ConstructEmptyDistributedSparseMatrix(V, exa%mat_dim)

    !! Get The Eigenvectors
    row_start = eigen_loop_start(1, exa%proc_rows, exa%rowid)
    row_end = eigen_loop_end(exa%mat_dim, exa%proc_rows, exa%rowid)
    col_start = eigen_loop_start(1, exa%proc_cols, exa%colid)
    col_end = eigen_loop_end(exa%mat_dim, exa%proc_cols, exa%colid)

    !! Convert to a 1D array for index ease.
    ALLOCATE(VD1(SIZE(VD,DIM=1)*SIZE(VD,DIM=2)))
    VD1 = PACK(VD, .TRUE.)

    CALL StartTimer("EigenExaFilter")
    CALL ConstructTripletList(triplet_v)
    ind = 1
    DO JJ = col_start, col_end
       jlookup = eigen_translate_l2g(JJ, exa%proc_cols, exa%colid)
       DO II = row_start, row_end
          IF (ABS(VD1(ind+II-1)) .GT. solver_parameters%threshold) THEN
             ilookup = eigen_translate_l2g(II, exa%proc_rows, exa%rowid)
             CALL SetTriplet(trip, jlookup, ilookup, VD1(ind+II-1))
             CALL AppendToTripletList(triplet_v, trip)
          END IF
       END DO
       ind = ind + exa%offset
    END DO
    CALL StopTimer("EigenExaFilter")

    CALL StartTimer("EigenFill")
    CALL FillFromTripletList(V, triplet_v)
    CALL StopTimer("EigenFill")

    !! Cleanup
    CALL DestructTripletList(triplet_w)
    CALL DestructTripletList(triplet_v)

    DEALLOCATE(VD1)

  END SUBROUTINE EigenToNT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The routine which calls the eigenexa driver.
  !! @param[in] A the matrix to decompose.
  !! @param[inout] V the eigenvectors.
  !! @param[inout] W the eigenvalues.
  !! @param[in] exa calculation parameters.
  SUBROUTINE Compute(A, V, W, exa)
    !! Parameters
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: A
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: V
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: W
    TYPE(ExaHelper_t), INTENT(IN) :: exa
    !! Local Variables
    INTEGER :: N, LDA, LDZ

    !! Setup EigenExa Parameters
    N = exa%mat_dim
    LDA = exa%local_rows
    LDZ = exa%local_rows

    !! Call
    CALL eigen_sx(N, N, A, LDA, W, V, LDZ, mode=exa%MODE)

  END SUBROUTINE Compute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocates and shuts down eigenexa
  !! @param[inout] AD the matrix.
  !! @param[inout] VD the eigenvectors.
  !! @param[inout] W the eigenvalues.
  SUBROUTINE CleanUp(AD, VD, WD)
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AD
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: VD
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WD

    IF(ALLOCATED(AD)) DEALLOCATE(AD)
    IF(ALLOCATED(VD)) DEALLOCATE(VD)
    IF(ALLOCATED(WD)) DEALLOCATE(WD)

    CALL eigen_free

  END SUBROUTINE CleanUp
END MODULE EigenExaModule
