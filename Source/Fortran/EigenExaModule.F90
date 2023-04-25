!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for calling eigenexa
MODULE EigenExaModule
#if EIGENEXA
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteHeader, WriteListElement
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & FillMatrixFromTripletList, GetMatrixTripletList, PrintMatrixInformation
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  USE TripletModule, ONLY : Triplet_r, Triplet_c, SetTriplet
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & AppendToTripletList, GetTripletAt, ConstructTripletList, &
       & DestructTripletList, RedistributeTripletLists
  USE eigen_libs_mod
  USE NTMPIModule
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
     !> The global rank
     INTEGER :: global_rank
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
     !> For block cyclic indexing.
     INTEGER :: offset
     !> Number of values to compute.
     INTEGER :: nvals
  END TYPE ExaHelper_t
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a matrix using EigenExa.
  SUBROUTINE EigenExa_s(A, eigenvalues, nvals, &
       & eigenvectors_in, solver_parameters_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The eigenvalues computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> The number of eigenvalues to compute.
    INTEGER, INTENT(IN) :: nvals
    !> The eigenvectors computed.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvectors_in
    !> The parameters for this solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Process Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    !! Write info about the solver
    IF (params%be_verbose) THEN
       CALL WriteHeader("Eigen Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="EigenExa")
       CALL WriteElement(key="NVALS", VALUE=nvals)
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("imamura2011development")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Select Based on Type
    IF (A%is_complex) THEN
       IF (PRESENT(eigenvectors_in)) THEN
          CALL EigenExa_c(A, eigenvalues, nvals, params, &
               & eigenvectors_in)
       ELSE
          CALL EigenExa_c(A, eigenvalues, nvals, params)
       END IF
    ELSE
       IF (PRESENT(eigenvectors_in)) THEN
          CALL EigenExa_r(A, eigenvalues, nvals, params, &
               & eigenvectors_in)
       ELSE
          CALL EigenExa_r(A, eigenvalues, nvals, params)
       END IF
    END IF

    !! Cleanup
    CALL DestructSolverParameters(params)

  END SUBROUTINE EigenExa_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a matrix using EigenExa (real).
  SUBROUTINE EigenExa_r(A, eigenvalues, nvals, params, eigenvectors_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The eigenvalues computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> The number of eigenvalues to compute.
    INTEGER, INTENT(IN) :: nvals
    !> The parameters for this solver.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !> The eigenvectors computed.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvectors_in
    !! Helper
    TYPE(ExaHelper_t) :: exa
    !! Dense Matrices
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: AD
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: VD
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: WD

#include "eigenexa_includes/EigenExa_s.F90"

  END SUBROUTINE EigenExa_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a matrix using EigenExa (complex).
  SUBROUTINE EigenExa_c(A, eigenvalues, nvals, params, eigenvectors_in)
    !> The matrix to decompose.
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The eigenvalues computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    !> The number of eigenvalues to compute.
    INTEGER, INTENT(IN) :: nvals
    !> The parameters for this solver.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !> The eigenvectors computed.
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: eigenvectors_in
    !! Helper
    TYPE(ExaHelper_t) :: exa
    !! Dense Matrices
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: AD
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: VD
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: WD

#define ISCOMPLEX
#include "eigenexa_includes/EigenExa_s.F90"
#undef ISCOMPLEX

  END SUBROUTINE EigenExa_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Setup the eigen exa data structures.
  SUBROUTINE InitializeEigenExa(A, nvals, eigenvectors, exa)
    !> The matrix we're working on.
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> Number of eigenvalues to compute.
    INTEGER, INTENT(IN) :: nvals
    !> Whether to compute eigenvectors.
    LOGICAL, INTENT(IN) :: eigenvectors
    !> Stores info about the calculation.
    TYPE(ExaHelper_t), INTENT(INOUT) :: exa
    !! Local Variables
    INTEGER :: ICTXT, INFO
    INTEGER, DIMENSION(9) :: DESCA
    INTEGER :: ierr

    !! Number of values to compute.
    exa%nvals = nvals

    !! Setup the MPI Communicator
    CALL MPI_Comm_dup(A%process_grid%global_comm, exa%comm, ierr)
    CALL MPI_Comm_rank(exa%comm, exa%global_rank, ierr)

    !! Build EigenExa Process Grid
    CALL eigen_init(exa%comm)
    CALL eigen_get_procs(exa%num_procs, exa%proc_rows, exa%proc_cols )
    CALL eigen_get_id(exa%procid, exa%rowid, exa%colid)

    !! Allocate Dense Matrices
    exa%mat_dim = A%actual_matrix_dimension
    CALL eigen_get_matdims(exa%mat_dim, exa%local_rows, exa%local_cols)

    !> Default blocking parameters
    exa%MB = 128
    exa%M = 48
    IF (eigenvectors) THEN
       exa%MODE = 'A'
    ELSE
       exa%MODE = 'N'
    END IF

    !! Blacs gives us the blocking info.
    ICTXT = eigen_get_blacs_context()
    CALL DESCINIT(DESCA, exa%mat_dim, exa%mat_dim, 1, 1, 0, 0, ICTXT, &
         & exa%local_rows, INFO )
    exa%offset = DESCA(9)

  END SUBROUTINE InitializeEigenExa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the distributed sparse matrix to a dense matrix in block-cyclic
  !> distribution (real).
  SUBROUTINE NTToEigen_r(A, AD, exa)
    !> The matrix to convert.
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The dense, block-cyclic version.
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: AD
    !> Info about the calculation.
    TYPE(ExaHelper_t), INTENT(INOUT) :: exa
    !! Local Variables
    TYPE(TripletList_r) :: triplet_a
    TYPE(TripletList_r), DIMENSION(:), ALLOCATABLE :: send_list
    TYPE(TripletList_r) :: recv_list
    TYPE(Triplet_r) :: trip, shifted_trip

#include "eigenexa_includes/NTToEigen.f90"

  END SUBROUTINE NTToEigen_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the distributed sparse matrix to a dense matrix in block-cyclic
  !> distribution (complex).
  SUBROUTINE NTToEigen_c(A, AD, exa)
    !> The matrix to convert.
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The dense, block-cyclic version.
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), INTENT(INOUT) :: AD
    !> Info about the calculation.
    TYPE(ExaHelper_t), INTENT(INOUT) :: exa
    !! Local Variables
    TYPE(TripletList_c) :: triplet_a
    TYPE(TripletList_c), DIMENSION(:), ALLOCATABLE :: send_list
    TYPE(TripletList_c) :: recv_list
    TYPE(Triplet_c) :: trip, shifted_trip

#include "eigenexa_includes/NTToEigen.f90"

  END SUBROUTINE NTToEigen_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the dense eigenvector matrix stored block-cyclicly back to
  !> a distributed sparse matrix (real).
  SUBROUTINE EigenToNT_r(VD, V, params, exa)
    !> The dense eigenvector matrix.
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: VD
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: V
    !> Parameters for thresholding small values.
    TYPE(SolverParameters_t) :: params
    !> Info about the calculation.
    TYPE(ExaHelper_t) :: exa
    !! Local Variables
    TYPE(TripletList_r) :: triplet_v
    TYPE(Triplet_r) :: trip
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: VD1

#include "eigenexa_includes/EigenToNT.f90"

  END SUBROUTINE EigenToNT_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the dense eigenvector matrix stored block-cyclicly back to
  !> a distributed sparse matrix (complex).
  SUBROUTINE EigenToNT_c(VD, V, params, exa)
    !> The dense eigenvector matrix.
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), INTENT(IN) :: VD
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: V
    !> Parameters for thresholding small values.
    TYPE(SolverParameters_t) :: params
    !> Info about the calculation.
    TYPE(ExaHelper_t) :: exa
    !! Local Variables
    TYPE(TripletList_c) :: triplet_v
    TYPE(Triplet_c) :: trip
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: VD1

#include "eigenexa_includes/EigenToNT.f90"

  END SUBROUTINE EigenToNT_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the dense eigenvalue matrix stored duplicated across processes.
  SUBROUTINE ExtractEigenvalues(WD, W, exa)
    !> The dense eigenvalue matrix.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: WD
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: W
    !> Info about the calculation.
    TYPE(ExaHelper_t) :: exa
    !! Local Variables
    TYPE(TripletList_r) :: triplet_w
    TYPE(Triplet_r) :: trip
    INTEGER :: wstart, wend, wsize
    INTEGER :: II

    !! Copy To Triplet List
    wsize = MAX(CEILING((1.0*exa%mat_dim)/exa%num_procs), 1)
    wstart = wsize*exa%global_rank + 1
    wend = MIN(wsize*(exa%global_rank+1), exa%mat_dim)

    CALL ConstructTripletList(triplet_w)
    DO II = wstart, wend
       IF (II .GT. exa%nvals) THEN
          EXIT
       END IF
       CALL SetTriplet(trip, II, II, WD(II))
       CALL AppendToTripletList(triplet_w, trip)
    END DO

    !! Go to global matrix
    CALL FillMatrixFromTripletList(W, triplet_w)

    !! Cleanup
    CALL DestructTripletList(triplet_w)

  END SUBROUTINE ExtractEigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The routine which calls the eigenexa driver.
  SUBROUTINE Compute_r(A, V, W, exa)
    !> The matrix to decompose.
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: A
    !> The eigenvectors.
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: V
    !> The eigenvalues.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: W
    !> Calculation parameters.
    TYPE(ExaHelper_t), INTENT(IN) :: exa

#include "eigenexa_includes/Compute.f90"

    !! Call
    CALL eigen_sx(N, exa%nvals, A, LDA, W, V, LDZ, mode=exa%MODE)

  END SUBROUTINE Compute_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The routine which calls the eigenexa driver.
  SUBROUTINE Compute_c(A, V, W, exa)
    !> The matrix to decompose.
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), INTENT(INOUT) :: A
    !> The eigenvectors.
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), INTENT(INOUT) :: V
    !> The eigenvalues.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: W
    !> Calculation parameters.
    TYPE(ExaHelper_t), INTENT(IN) :: exa

#include "eigenexa_includes/Compute.f90"

    !! Call
    CALL eigen_h(N, exa%nvals, A, LDA, W, V, LDZ, mode=exa%MODE)

  END SUBROUTINE Compute_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocates and shuts down eigenexa (real)
  SUBROUTINE CleanUp_r(AD, VD, WD)
    !> The matrix._r
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AD
    !> The eigenvectors.
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: VD
    !> The eigenvalues.
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WD

#include "eigenexa_includes/Cleanup.f90"

  END SUBROUTINE CleanUp_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocates and shuts down eigenexa (complex)
  SUBROUTINE CleanUp_c(AD, VD, WD)
    !> The matrix.
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AD
    !> The eigenvectors.
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: VD
    !> The eigenvalues.
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WD

#include "eigenexa_includes/Cleanup.f90"

  END SUBROUTINE CleanUp_c
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenExaModule
