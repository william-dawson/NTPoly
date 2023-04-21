!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Density Matrix Using the Fermi Operator Expansion
MODULE FermiOperatorModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE EigenSolversModule, ONLY : EigenDecomposition
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteHeader, &
       & EnterSubLog, ExitSubLog
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & FillMatrixFromTripletList, GetMatrixTripletList, &
       & TransposeMatrix, ConjugateMatrix, DestructMatrix
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE SolverParametersModule, ONLY : SolverParameters_t, &
       & PrintParameters, DestructSolverParameters, CopySolverParameters, &
       & ConstructSolverParameters
  USE TripletListModule, ONLY : TripletList_r, DestructTripletList
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeDenseFOE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix using a dense routine.
  SUBROUTINE ComputeDenseFOE(H, ISQ, trace, K, inv_temp_in, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The inverse temperature for smearing (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: inv_temp_in
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    REAL(NTREAL) :: inv_temp
    LOGICAL :: do_smearing
    !! Local Variables
    TYPE(Matrix_ps) :: ISQT, WH
    TYPE(Matrix_ps) :: WD
    TYPE(Matrix_ps) :: vecs, vecsT, vals, Temp
    TYPE(MatrixMemoryPool_p) :: pool
    TYPE(TripletList_r) :: tlist
    REAL(NTREAL) :: chemical_potential, energy_value
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: eigs, occ
    REAL(NTREAL) :: sval, sv, occ_temp
    REAL(NTREAL) :: left, right, homo, lumo
    INTEGER :: num_eigs
    INTEGER :: II, JJ
    INTEGER :: ierr

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF
    IF (PRESENT(inv_temp_in)) THEN
       inv_temp = inv_temp_in
       do_smearing = .TRUE.
    ELSE
       do_smearing = .FALSE.
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       IF (do_smearing) THEN
          CALL WriteElement(key="Method", VALUE="Dense FOE")
          CALL WriteElement(key="Inverse Temperature", VALUE=inv_temp)
       ELSE
          CALL WriteElement(key="Method", VALUE="Dense Step Function")
       END IF
       CALL PrintParameters(params)
    END IF

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL MatrixMultiply(ISQ, H, Temp, &
         & threshold_in=params%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(Temp, ISQT, WH, &
         & threshold_in=params%threshold, memory_pool_in=pool)

    !! Perform the eigendecomposition
    CALL EigenDecomposition(WH, vals, &
         & eigenvectors_in=vecs, solver_parameters_in=params)

    !! Gather the eigenvalues on to every process
    CALL GetMatrixTripletList(vals, tlist)
    num_eigs = H%actual_matrix_dimension
    ALLOCATE(eigs(num_eigs))
    eigs = 0
    DO II = 1, tlist%CurrentSize
       eigs(tlist%DATA(II)%index_column) = tlist%DATA(II)%point_value
    END DO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, eigs, num_eigs, MPINTREAL, &
         & MPI_SUM, H%process_grid%within_slice_comm, ierr)

    !! Compute MU By Bisection
    IF (do_smearing) THEN
       ALLOCATE(occ(num_eigs))
       left = MINVAL(eigs)
       right = MAXVAL(eigs)
       DO JJ = 1, 10*params%max_iterations 
          chemical_potential = left + (right - left)/2 
          DO II = 1, num_eigs
             sval = eigs(II) - chemical_potential
             occ(II) = 0.5_NTREAL * (1.0_NTREAL - ERF(inv_temp * sval))
          END DO
          sv = SUM(occ)
          IF (ABS(trace - sv) .LT. 1E-8_NTREAL) THEN
             EXIT
          ELSE IF (SV > trace) THEN
             right = chemical_potential
          ELSE
             left = chemical_potential
          END IF
       END DO
    ELSE
       JJ = 1
       homo = eigs(FLOOR(trace))
       lumo = eigs(FLOOR(trace) + 1)
       occ_temp = FLOOR(TRACE) + 1 - trace
       chemical_potential = homo + occ_temp * 0.5_NTREAL * (lumo - homo)
    END IF

    !! Write out result of chemical potential search
    IF (params%be_verbose) THEN
       CALL WriteHeader("Chemical Potential Search")
       CALL EnterSubLog
       CALL WriteElement(key="Potential", VALUE=chemical_potential)
       CALL WriteElement(key="Iterations", VALUE=JJ)
       CALL ExitSubLog
    END IF

    !! Map
    energy_value = 0.0_NTREAL
    DO II = 1, tlist%CurrentSize
       IF (.NOT. do_smearing) THEN
          IF (tlist%DATA(II)%index_column .LE. FLOOR(trace)) THEN
             energy_value = energy_value + tlist%DATA(II)%point_value
             tlist%DATA(II)%point_value = 1.0_NTREAL
          ELSE IF (tlist%DATA(II)%index_column .EQ. CEILING(trace)) THEN
             occ_temp = CEILING(trace) - trace
             energy_value = energy_value + &
                  & occ_temp * tlist%DATA(II)%point_value
             tlist%DATA(II)%point_value = occ_temp
          ELSE
             tlist%DATA(II)%point_value = 0.0_NTREAL
          ENDIF
       ELSE
          sval = tlist%DATA(II)%point_value - chemical_potential
          occ_temp = 0.5_NTREAL * (1.0_NTREAL - ERF(inv_temp * sval))
          energy_value = energy_value + occ_temp * tlist%DATA(II)%point_value
          tlist%DATA(II)%point_value = occ_temp
       END IF
    END DO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, energy_value, 1, MPINTREAL, MPI_SUM, &
         & H%process_grid%within_slice_comm, ierr)

    !! Fill
    CALL ConstructEmptyMatrix(vals, H)
    CALL FillMatrixFromTripletList(vals, tlist, preduplicated_in=.TRUE.)

    !! Multiply Back Together
    CALL MatrixMultiply(vecs, vals, temp, threshold_in=params%threshold)
    CALL TransposeMatrix(vecs, vecsT)
    CALL ConjugateMatrix(vecsT)
    CALL MatrixMultiply(temp, vecsT, WD, &
         & threshold_in=params%threshold)

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(ISQT, WD, Temp, &
         & threshold_in=params%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(Temp, ISQ, K, &
         & threshold_in=params%threshold, memory_pool_in=pool)

    !! Optional out variables.
    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = 2.0_NTREAL * energy_value
    END IF
    IF (PRESENT(chemical_potential_out)) THEN
       chemical_potential_out = chemical_potential
    END IF

    !! Cleanup
    CALL DestructMatrix(WH)
    CALL DestructMatrix(WD)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(vecs)
    CALL DestructMatrix(vecst)
    CALL DestructMatrix(vals)
    CALL DestructMatrix(temp)
    CALL DestructTripletList(tlist)
    CALL DestructMatrixMemoryPool(pool)
    IF (ALLOCATED(occ)) THEN
       DEALLOCATE(occ)
    END IF
    IF (ALLOCATED(eigs)) THEN
       DEALLOCATE(eigs)
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
  END SUBROUTINE ComputeDenseFOE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorModule
