!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Density Matrix Using the Fermi Operator Expansion
MODULE FermiOperatorExpansionModule
  USE ChebyshevSolversModule, ONLY : ChebyshevPolynomial_t, &
       & ConstructPolynomial, DestructPolynomial, SetCoefficient, &
       & Compute
  USE DataTypesModule, ONLY : NTREAL
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteHeader, EnterSubLog, &
       & ExitSubLog
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, FillMatrixIdentity, &
       & PrintMatrixInformation, ConstructEmptyMatrix, DestructMatrix, &
       & CopyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeFOE
  PUBLIC :: FOEEigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> PI
  REAL(NTREAL), PARAMETER, PRIVATE :: PI =  4 * ATAN (1.0_16)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Density Matrix using the Fermi Operator Expansion.
  SUBROUTINE ComputeFOE(Hamiltonian, InverseSquareRoot, nel, Density, degree, &
       & chemical_potential_in, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN)  :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> Degree to expand to.
    INTEGER, INTENT(IN) :: degree
    !> If the chemical potential is already known, you can pass it in.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: chemical_potential_in
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: TempMat
    !! Local Variables
    REAL(NTREAL), DIMENSION(degree+1) :: moments, damping
    TYPE(ChebyshevPolynomial_t) :: cheb
    INTEGER :: II

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    !! Compute The Working Hamiltonian
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="FOE")
       CALL WriteElement(key="Degree", int_value_in=degree)
       CALL PrintParameters(solver_parameters)
    END IF

    !! If we know the chemical potential, do the standard algorithm.
    IF (PRESENT(chemical_potential_in)) THEN
       CALL ComputeFixed(WorkingHamiltonian, Density, degree, &
            & chemical_potential_in, solver_parameters)
    ELSE

    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(InverseSquareRoot,Density,TempMat, &
         & threshold_in=solver_parameters%threshold)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold)

    !! Cleanup
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(WorkingHamiltonian)

  END SUBROUTINE ComputeFOE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Estimate the eigenvalues of a matrix using the Fermi Operator Expansion.
  SUBROUTINE FOEEigenvalues(InputMat, eigenvalues, degree, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The eigenvalues of the matrix, preallocated with a length of the number
    !> of eigenvalues you want to compute.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: eigenvalues
    !> Degree to expand to.
    INTEGER, INTENT(IN) :: degree
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
  END SUBROUTINE FOEEigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Implementation with a fixed chemical potential.
  SUBROUTINE ComputeFixed(Hamiltonian, Density, degree, &
       & chemical_potential, solver_parameters)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: Hamiltonian
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> Degree to expand to.
    INTEGER, INTENT(IN) :: degree
    !> If the chemical potential is already known, you can pass it in.
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
    !! Local Variables
    REAL(NTREAL), DIMENSION(degree+1) :: moments, damping
    TYPE(ChebyshevPolynomial_t) :: cheb
    INTEGER :: II

    CALL ComputeMoments(moments, chemical_potential)
    CALL ComputeDamping(damping)

    CALL ConstructPolynomial(cheb, degree+1)
    DO II = 1, degree+1
       CALL SetCoefficient(cheb, II, moments(II)*damping(II))
    END DO
    CALL Compute(Hamiltonian, Density, cheb, solver_parameters)
    CALL DestructPolynomial(cheb)
    !! Otherwise, we have to do a search.


  END SUBROUTINE ComputeFixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the chebyshev moments to expand the fermi operator.
  !> Equation 1 in \cite{jay1999electronic}.
  SUBROUTINE ComputeMoments(moments, chemical_potential)
    !> The moments to compute.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: moments
    !> The chemical potential
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    !! Local Variables
    INTEGER :: II
    INTEGER :: M
    REAL(NTREAL) :: denom, numer

    !! Get the dimensions
    M = SIZE(moments) - 1

    !! Actual Loop
    moments(1)  = 2.0 * (1.0 - ACOS(chemical_potential)/PI)
    DO II = 2, M + 1
       denom = (II-1)*PI
       numer = 2.0*SIN((II-1)*ACOS(chemical_potential))
       moments(II) = -1.0 * numer/denom
    END DO
  END SUBROUTINE ComputeMoments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the damping factor to remove the gibbs oscillations.
  !> Equation 5 in \cite{jay1999electronic}.
  SUBROUTINE ComputeDamping(damping)
    !> The damping coefficients to compute.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: damping
    !! Local Variables
    REAL(NTREAL) :: left, right, denom
    INTEGER :: II
    INTEGER :: M
    INTEGER :: Mplus2

    !! Get the dimensions
    M = SIZE(damping) - 1
    Mplus2 = M + 2

    !! Actual Loop
    DO II = 1, M + 1
       left = (1 - (II-1)/(Mplus2))*SIN(PI/Mplus2)*COS(((II-1)*PI)/Mplus2)
       right = (1.0/Mplus2)*COS(PI/Mplus2)*SIN(((II-1)*PI)/Mplus2)
       denom = SIN(PI/Mplus2)
       damping(II) = (left + right)/denom
    END DO
  END SUBROUTINE ComputeDamping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the threshold estimates. We do this by taking an arbitrary
  !> value for the chemical potential, and then smoothing the threshold
  !> values to make sure they apply to other potential values.
  SUBROUTINE ComputeThresholds(estimates, damping, threshold)
    !> The estimates we're computing.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: estimates
    !> The damping factors.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: damping
    !> The threshold parameter.
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local variables
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: moments
    INTEGER :: bstart, bend, bsize
    INTEGER :: II, degree

    degree = SIZE(estimates)

    !! Compute the threshold values using an initial guess of moments
    ALLOCATE(moments(degree))
    CALL ComputeMoments(moments, REAL(0.25, NTREAL))
    estimates = ABS(threshold * moments*damping)

    !! Using the max value in a range, we can get a good convervative bounds.
    bsize = 16
    bend = 1

    !! Main loop
    DO II = 0, (degree)/bsize - 1
       bstart = II * bsize + 1
       bend = (II+1) * bsize
       estimates(bstart:bend) = MAXVAL(estimates(bstart:bend))
    END DO
    !! Get the final values
    estimates(bend:) = MAXVAL(estimates(bend:))

    !! Cleanup
    DEALLOCATE(moments)
  END SUBROUTINE ComputeThresholds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorExpansionModule
