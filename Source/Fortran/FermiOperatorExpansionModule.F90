!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Density Matrix Using the Fermi Operator Expansion
MODULE FermiOperatorExpansionModule
  USE ChebyshevSolversModule, ONLY : ChebyshevPolynomial_t, &
       & ConstructPolynomial, DestructPolynomial, SetCoefficient, &
       & FactorizedCompute
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteHeader, EnterSubLog, &
       & ExitSubLog
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, &
       & ScaleMatrix, MatrixTrace, DotMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, FillMatrixIdentity, &
       & PrintMatrixInformation, ConstructEmptyMatrix, DestructMatrix, &
       & CopyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeFOE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> PI
  REAL(NTREAL), PARAMETER, PRIVATE :: PI =  4 * ATAN (1.0_NTREAL)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute The Density Matrix using the Fermi Operator Expansion.
  SUBROUTINE ComputeFOE(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN)  :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The computed chemical potential.
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: TempMat
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: energy_value, chemical_potential

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

    !! Scale the matrix.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)
    CALL ScaleMatrix(WorkingHamiltonian, REAL(1.0,NTREAL)/(e_max-e_min))

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="FOE")
       CALL WriteElement(key="Degree", value=solver_parameters%max_iterations)
       CALL PrintParameters(solver_parameters)
    END IF

    !! If we know the chemical potential, do the standard algorithm.
    CALL ComputeSearch(WorkingHamiltonian, Density, nel/2, &
         & solver_parameters%max_iterations, &
         & solver_parameters, energy_value, chemical_potential)
    IF (PRESENT(chemical_potential_out)) THEN
       chemical_potential_out = chemical_potential * (e_max - e_min)
    END IF
    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
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
  !> Implementation that has to search for the chemical potential
  SUBROUTINE ComputeSearch(Hamiltonian, Density, target_trace, degree, &
       & solver_parameters, energy_value, chemical_potential)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: Hamiltonian
    !> OutputMat = poly(InputMat)
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> Target trace of the matrix.
    INTEGER, INTENT(IN) :: target_trace
    !> Degree to expand to.
    INTEGER, INTENT(IN) :: degree
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
    !> The energy of the system.
    REAL(NTREAL), INTENT(OUT) :: energy_value
    !> The computed chemical potential.
    REAL(NTREAL), INTENT(OUT) :: chemical_potential
    !! Local Variables
    TYPE(Matrix_ps), DIMENSION(:), ALLOCATABLE :: PowerArray
    REAL(NTREAL), DIMENSION(degree+1) :: damping, thresholds
    REAL(NTREAL) :: cp_left, cp_right
    REAL(NTREAL) :: trace_value
    INTEGER :: II
    LOGICAL :: finished

    !! Make the coefficients
    CALL ComputeDamping(damping)
    CALL ComputeThresholds(thresholds, damping, solver_parameters%threshold)

    !! Core computation of the power arrays
    ALLOCATE(PowerArray(degree))
    CALL ComputePowers(Hamiltonian, PowerArray, thresholds, solver_parameters)

    !! Search
    cp_left = -1.0
    cp_right = 1.0
    finished = .FALSE.
    DO WHILE(.NOT. finished)
       chemical_potential = cp_left + (cp_right - cp_left)/2.0
       CALL SumFOE(PowerArray, damping, chemical_potential, degree, Density)
       CALL MatrixTrace(Density, trace_value)

       IF (ABS(trace_value - target_trace) .LT. &
            & solver_parameters%converge_diff) THEN
          finished = .TRUE.
       ELSE IF (trace_value .GT. target_trace) THEN
          cp_right = chemical_potential
       ELSE
          cp_left = chemical_potential
       END IF
    END DO
    CALL DotMatrix(Hamiltonian, Density, energy_value)

    !! Undo load balancing
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(Density, Density, &
            & solver_parameters%BalancePermutation)
    END IF

    !! Cleanup
    DO II = 1, degree
       CALL DestructMatrix(PowerArray(II))
    END DO
    DEALLOCATE(PowerArray)

  END SUBROUTINE ComputeSearch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the powers of the Chebysehv polynomial
  !> note that this routine won't reverse the permutation.
  SUBROUTINE ComputePowers(InputMat, PowerArray, thresholds, solver_parameters)
    !> The starting matrix
    TYPE(Matrix_ps), INTENT(IN) :: InputMat
    !> An array of Chebysehv polynomials to compute
    TYPE(Matrix_ps), DIMENSION(:), INTENT(INOUT) :: PowerArray
    !> Thresholds for each power
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: thresholds
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: solver_parameters
    !! Local Variables
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: BalancedInput
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II, degree

    !! Basic Setup
    degree = SIZE(PowerArray)
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL FillMatrixIdentity(Identity)
    CALL CopyMatrix(InputMat,BalancedInput)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(BalancedInput, BalancedInput, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! First matrices
    IF (degree .GT. 0) THEN
       CALL CopyMatrix(Identity, PowerArray(1))
    END IF
    IF (degree .GT. 1) THEN
       CALL CopyMatrix(BalancedInput, PowerArray(2))
    END IF

    !! Full Loop
    DO II = 3, degree
       CALL MatrixMultiply(BalancedInput, PowerArray(II-1), PowerArray(II), &
            & alpha_in=REAL(2.0,NTREAL), threshold_in=thresholds(II), &
            & memory_pool_in=pool)
       CALL IncrementMatrix(PowerArray(II-2), PowerArray(II), REAL(-1.0,NTREAL))
    END DO

    !! Cleanup
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(BalancedInput)
  END SUBROUTINE ComputePowers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the matrices for a given chemical potential
  SUBROUTINE SumFOE(PowerArray, damping, chemical_potential, degree, Density)
    !> The Chebyshev polynomials to sum up.
    TYPE(Matrix_ps), DIMENSION(:), INTENT(INOUT) :: PowerArray
    !> The damping factor.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: damping
    !> The chemical potential we are guessing.
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    !> The degree of the polynomial.
    INTEGER, INTENT(IN) :: degree
    !> The density value computed.
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !! Local variables
    INTEGER :: II
    REAL(NTREAL), DIMENSION(degree) :: moments

    !! Compute coefficients
    CALL ComputeMoments(moments,chemical_potential)

    !! Sum Up
    CALL CopyMatrix(PowerArray(1), Density)
    CALL ScaleMatrix(Density, moments(1)*damping(1))
    DO II = 2, degree
       CALL IncrementMatrix(PowerArray(II), Density, &
            & alpha_in=moments(II)*damping(II))
    END DO

    !! Final increment
    CALL IncrementMatrix(PowerArray(1), Density, alpha_in=-0.5*moments(1))
  END SUBROUTINE SumFOE
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
       left = (1 - (II-1.0)/(Mplus2))*SIN(PI/Mplus2)*COS(((II-1.0)*PI)/Mplus2)
       right = (1.0/Mplus2)*COS(PI/Mplus2)*SIN(((II-1.0)*PI)/Mplus2)
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
    estimates = threshold

    !! Cleanup
    DEALLOCATE(moments)
  END SUBROUTINE ComputeThresholds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorExpansionModule
