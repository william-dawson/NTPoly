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
       & PrintParameters, DestructSolverParameters
  USE TripletListModule, ONLY : TripletList_r, DestructTripletList
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeDenseFOE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix using a dense routine.
  SUBROUTINE ComputeDenseFOE(H, ISQ, nel, K, inv_temp_in, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
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
    REAL(NTREAL) :: homo, lumo
    REAL(NTREAL) :: sval, occ
    INTEGER :: II
    INTEGER :: ierr

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       params = solver_parameters_in
    ELSE
       params = SolverParameters_t()
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
          CALL WriteElement(key="InverseTemperature", VALUE=inv_temp)
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
    CALL EigenDecomposition(WH, vecs, vals, params)

    !! Convert to a triplet list and get homo/lumo + energy.
    CALL GetMatrixTripletList(vals, tlist)
    homo = 0.0_NTREAL
    lumo = 0.0_NTREAL
    energy_value = 0.0_NTREAL
    DO II = 1, tlist%CurrentSize
       IF (tlist%DATA(II)%index_column .EQ. INT(nel/2)) THEN
          homo = tlist%DATA(II)%point_value
       ELSE IF (tlist%DATA(II)%index_column .EQ. INT(nel/2) + 1) THEN
          lumo = tlist%DATA(II)%point_value
       END IF
    END DO

    !! Compute MU
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, homo, 1, MPINTREAL, MPI_SUM, &
         & H%process_grid%within_slice_comm, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, lumo, 1, MPINTREAL, MPI_SUM, &
         & H%process_grid%within_slice_comm, ierr)
    chemical_potential = homo + 0.5_NTREAL * (lumo - homo)

    !! Map
    DO II = 1, tlist%CurrentSize
       IF (.NOT. do_smearing) THEN
          IF (tlist%DATA(II)%index_column .LE. INT(nel/2)) THEN
             energy_value = energy_value + tlist%DATA(II)%point_value
             tlist%DATA(II)%point_value = 1.0_NTREAL
          ELSE
             tlist%DATA(II)%point_value = 0.0_NTREAL
          ENDIF
       ELSE
          sval = tlist%DATA(II)%point_value - chemical_potential
          occ = 1.0 / (1 + EXP(inv_temp * sval))
          energy_value = energy_value + occ * tlist%DATA(II)%point_value
          tlist%DATA(II)%point_value = occ
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

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
  END SUBROUTINE ComputeDenseFOE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorModule
