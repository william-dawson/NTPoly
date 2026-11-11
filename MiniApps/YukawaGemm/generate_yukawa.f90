!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Generate a synthetic sparse matrix from a Yukawa potential on a simple cubic
!> lattice, writing MatrixMarket format. Self-contained Fortran (no scipy,
!> no MPI), so it compiles and runs anywhere.
!>
!> V_ij = exp(-r_ij / lambda) / r_ij  for r_ij <= Rcut, else 0
!> diagonal = diag (default 1.0)
!>
!> Neighbor search uses a lattice-offset stencil — O(N * stencil^3),
!> no KD-tree needed for a regular lattice.
PROGRAM GenerateYukawa
  IMPLICIT NONE
  INTEGER :: side, seed
  REAL(8) :: a, lambda, rcut, jitter, diag
  CHARACTER(LEN=1024) :: outfile
  INTEGER :: nmax, ox, oy, oz, i, j, k, ni, nj, nk, atom_i, atom_j
  INTEGER(8) :: nnz
  REAL(8) :: dx, dy, dz, dist, val
  REAL(8), ALLOCATABLE :: pos(:,:)
  INTEGER :: unit_handle, nargs, ia
  CHARACTER(LEN=32) :: arg, argval
  INTEGER :: seed_size
  INTEGER, ALLOCATABLE :: seed_arr(:)

  !! Defaults matching the Python generator presets
  side = 32
  a = 2.7d0
  lambda = 1.0d0
  rcut = 6.0d0
  jitter = 0.05d0
  diag = 1.0d0
  seed = 42
  outfile = "yukawa.mtx"

  !! Parse arguments
  ia = 1
  nargs = COMMAND_ARGUMENT_COUNT()
  DO WHILE (ia .LE. nargs)
    CALL GET_COMMAND_ARGUMENT(ia, arg)
    CALL GET_COMMAND_ARGUMENT(ia+1, argval)
    SELECT CASE (TRIM(arg))
    CASE ('--side');      READ(argval,*) side
    CASE ('--lattice');   READ(argval,*) a
    CASE ('--screening'); READ(argval,*) lambda
    CASE ('--rcut');      READ(argval,*) rcut
    CASE ('--jitter');    READ(argval,*) jitter
    CASE ('--diag');      READ(argval,*) diag
    CASE ('--seed');      READ(argval,*) seed
    CASE ('--out');       outfile = TRIM(argval)
    END SELECT
    ia = ia + 2
  END DO

  !! Seed the RNG
  CALL RANDOM_SEED(SIZE=seed_size)
  ALLOCATE(seed_arr(seed_size))
  seed_arr = seed
  CALL RANDOM_SEED(PUT=seed_arr)
  DEALLOCATE(seed_arr)

  !! Build positions
  ALLOCATE(pos(3, side*side*side))
  atom_i = 0
  DO i = 0, side-1
    DO j = 0, side-1
      DO k = 0, side-1
        atom_i = atom_i + 1
        pos(1, atom_i) = REAL(i, 8) * a
        pos(2, atom_i) = REAL(j, 8) * a
        pos(3, atom_i) = REAL(k, 8) * a
        IF (jitter .GT. 0.0d0) THEN
          pos(1, atom_i) = pos(1, atom_i) + gauss() * jitter * a
          pos(2, atom_i) = pos(2, atom_i) + gauss() * jitter * a
          pos(3, atom_i) = pos(3, atom_i) + gauss() * jitter * a
        END IF
      END DO
    END DO
  END DO

  nmax = INT(rcut / a) + 1

  !! First pass: count nnz
  nnz = 0
  DO atom_i = 1, side*side*side
    CALL ijk_from_atom(atom_i, side, i, j, k)
    DO ox = -nmax, nmax
      DO oy = -nmax, nmax
        DO oz = -nmax, nmax
          IF (ox .EQ. 0 .AND. oy .EQ. 0 .AND. oz .EQ. 0) CYCLE
          ni = i + ox; nj = j + oy; nk = k + oz
          IF (ni < 0 .OR. ni >= side) CYCLE
          IF (nj < 0 .OR. nj >= side) CYCLE
          IF (nk < 0 .OR. nk >= side) CYCLE
          atom_j = atom_from_ijk(ni, nj, nk, side)
          IF (atom_j <= atom_i) CYCLE
          dx = pos(1, atom_i) - pos(1, atom_j)
          dy = pos(2, atom_i) - pos(2, atom_j)
          dz = pos(3, atom_i) - pos(3, atom_j)
          dist = SQRT(dx*dx + dy*dy + dz*dz)
          IF (dist > rcut .OR. dist < 1.0d-10) CYCLE
          nnz = nnz + 2
        END DO
      END DO
    END DO
  END DO
  nnz = nnz + side*side*side  ! diagonal

  WRITE(*,'(A,I0)') "atoms      : ", side*side*side
  WRITE(*,'(A,I0)') "nnz        : ", nnz
  WRITE(*,'(A,F0.1)') "nnz/row    : ", REAL(nnz, 8) / REAL(side*side*side, 8)

  !! Write MatrixMarket file (two-pass: header known from count)
  OPEN(NEWUNIT=unit_handle, FILE=TRIM(outfile), STATUS='REPLACE', ACTION='WRITE')
  WRITE(unit_handle,'(A)') "%%MatrixMarket matrix coordinate real general"
  WRITE(unit_handle,'(I0,1X,I0,1X,I0)') side*side*side, side*side*side, nnz

  !! Diagonal entries
  DO atom_i = 1, side*side*side
    WRITE(unit_handle,'(I0,1X,I0,1X,ES24.16E3)') atom_i, atom_i, diag
  END DO

  !! Off-diagonal entries (both directions for each pair)
  DO atom_i = 1, side*side*side
    CALL ijk_from_atom(atom_i, side, i, j, k)
    DO ox = -nmax, nmax
      DO oy = -nmax, nmax
        DO oz = -nmax, nmax
          IF (ox .EQ. 0 .AND. oy .EQ. 0 .AND. oz .EQ. 0) CYCLE
          ni = i + ox; nj = j + oy; nk = k + oz
          IF (ni < 0 .OR. ni >= side) CYCLE
          IF (nj < 0 .OR. nj >= side) CYCLE
          IF (nk < 0 .OR. nk >= side) CYCLE
          atom_j = atom_from_ijk(ni, nj, nk, side)
          IF (atom_j <= atom_i) CYCLE
          dx = pos(1, atom_i) - pos(1, atom_j)
          dy = pos(2, atom_i) - pos(2, atom_j)
          dz = pos(3, atom_i) - pos(3, atom_j)
          dist = SQRT(dx*dx + dy*dy + dz*dz)
          IF (dist > rcut .OR. dist < 1.0d-10) CYCLE
          val = EXP(-dist / lambda) / dist
          WRITE(unit_handle,'(I0,1X,I0,1X,ES24.16E3)') atom_i, atom_j, val
          WRITE(unit_handle,'(I0,1X,I0,1X,ES24.16E3)') atom_j, atom_i, val
        END DO
      END DO
    END DO
  END DO
  CLOSE(unit_handle)

  WRITE(*,'(A,A)') "wrote      : ", TRIM(outfile)
  DEALLOCATE(pos)

CONTAINS
  INTEGER FUNCTION atom_from_ijk(i, j, k, s)
    INTEGER, INTENT(IN) :: i, j, k, s
    atom_from_ijk = i * s * s + j * s + k + 1
  END FUNCTION
  SUBROUTINE ijk_from_atom(atom, s, i, j, k)
    INTEGER, INTENT(IN) :: atom, s
    INTEGER, INTENT(OUT) :: i, j, k
    INTEGER :: rem
    i = (atom - 1) / (s * s)
    rem = MOD(atom - 1, s * s)
    j = rem / s
    k = MOD(rem, s)
  END SUBROUTINE
  FUNCTION gauss() RESULT(g)
    REAL(8) :: g, u1, u2
    CALL RANDOM_NUMBER(u1)
    CALL RANDOM_NUMBER(u2)
    g = SQRT(-2.0d0 * LOG(MAX(u1, 1.0d-15))) * COS(2.0d0 * 3.14159265358979d0 * u2)
  END FUNCTION
END PROGRAM GenerateYukawa
