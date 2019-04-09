'''@package test_chemistry
A test suite for chemistry focused solvers.
'''
from helpers import THRESHOLD, EXTRAPTHRESHOLD
from helpers import result_file
from helpers import scratch_dir
import unittest
import NTPolySwig as nt
from numpy import diag, sqrt
from scipy.sparse.linalg import norm
from scipy.io import mmread, mmwrite
from scipy.linalg import eigh, funm
from scipy.sparse import csr_matrix, rand
from os import environ
from os.path import join
from mpi4py import MPI
# MPI global communicator
comm = MPI.COMM_WORLD


class TestChemistry:
    '''A test class for the chemistry solvers.'''
    # Matrix to compare to
    CheckMat = 0
    # Rank of the current process
    my_rank = 0

    def create_matrices(self):
        '''
        Create the test matrix with the following parameters.
        '''
        fock = rand(self.mat_dim, self.mat_dim, density=1.0)
        if self.is_complex:
            fock += 1j * rand(self.mat_dim, self.mat_dim, density=1.0)
            fock = fock + fock.H
        else:
            fock = fock + fock.T
        overlap = rand(self.mat_dim, self.mat_dim, density=1.0)
        overlap = overlap.T.dot(overlap)

        # Make sure the overlap is well conditioned.
        w, v = eigh(overlap.todense())
        w += 0.2
        overlap = csr_matrix(v.T.dot(diag(w).dot(v)))

        isq = funm(overlap.todense(), lambda x: 1.0 / sqrt(x))
        wfock = isq.dot(fock.todense()).dot(isq)

        # Add a gap
        w, v = eigh(wfock)
        gap = (w[-1] - w[0]) / 2.0
        w[self.nel:] += gap
        if self.is_complex:
            wfock = v.conj().T.dot(diag(w).dot(v))
        else:
            wfock = v.T.dot(diag(w).dot(v))

        # Compute the density
        w[:int(self.nel / 2)] = 2.0
        w[int(self.nel / 2):] = 0.0
        if self.is_complex:
            density = isq.dot(v.dot(diag(w).dot(v.conj().T))).dot(isq)
        else:
            density = isq.dot(v.dot(diag(w).dot(v.T))).dot(isq)

        self.write_matrix(fock, self.hamiltonian)
        self.write_matrix(overlap, self.overlap)
        self.write_matrix(density, self.density)

    def write_matrix(self, mat, file_name):
        if self.my_rank == 0:
            mmwrite(file_name, csr_matrix(mat))
        comm.barrier()

    @classmethod
    def setUpClass(self):
        '''Set up all of the tests.'''
        rows = int(environ['PROCESS_ROWS'])
        columns = int(environ['PROCESS_COLUMNS'])
        slices = int(environ['PROCESS_SLICES'])
        nt.ConstructGlobalProcessGrid(rows, columns, slices, True)

    @classmethod
    def tearDownClass(self):
        '''Cleanup this test'''
        nt.DestructGlobalProcessGrid()

    def setUp(self):
        '''Set up an individual test.'''
        self.my_rank = comm.Get_rank()
        self.solver_parameters = nt.SolverParameters()
        self.solver_parameters.SetVerbosity(True)
        self.geomh1 = environ["GEOMH1"]
        self.geomo1 = environ["GEOMO1"]
        self.geomo2 = environ["GEOMO2"]
        self.geomd2 = environ["GEOMD2"]
        self.realio = environ["REALIO"]
        self.nel = 10

        self.hamiltonian = join(scratch_dir, "rf.mtx")
        self.overlap = join(scratch_dir, "rs.mtx")
        self.density = join(scratch_dir, "rd.mtx")
        self.mat_dim = 7

    def check_full(self):
        '''Compare two computed matrices.'''
        normval = 0
        if (self.my_rank == 0):
            ResultMat = 2.0 * mmread(result_file)
            self.CheckMat = mmread(self.density)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    def check_full_extrap(self):
        '''Compare two computed matrices.'''
        normval = 0
        if (self.my_rank == 0):
            ResultMat = 2.0 * mmread(result_file)
            self.CheckMat = mmread(self.geomd2)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, EXTRAPTHRESHOLD)

    def compute_cp(self):
        '''Compute the chemical potential'''
        fock_matrix = mmread(self.hamiltonian)
        overlap_matrix = mmread(self.overlap)
        eig_vals = eigh(a=fock_matrix.todense(), b=overlap_matrix.todense(),
                        eigvals_only=True)
        homo = int(self.nel / 2) - 1
        lumo = homo + 1
        return eig_vals[homo] + (eig_vals[lumo] - eig_vals[homo]) / 2.0

    def check_cp(self, computed):
        '''Compare two computed chemical potentials.'''
        fock_matrix = mmread(self.hamiltonian)
        overlap_matrix = mmread(self.overlap)
        eig_vals = eigh(a=fock_matrix.todense(), b=overlap_matrix.todense(),
                        eigvals_only=True)
        homo = int(self.nel / 2) - 1
        lumo = homo + 1

        if computed > eig_vals[homo] and computed < eig_vals[lumo]:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def basic_solver(self, SRoutine, cpcheck=True):
        '''Test various kinds of density matrix solvers.'''
        self.create_matrices()
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        overlap_matrix = nt.Matrix_ps(self.overlap)
        inverse_sqrt_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        density_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        energy_value, chemical_potential = SRoutine(fock_matrix,
                                                    inverse_sqrt_matrix,
                                                    self.nel, density_matrix,
                                                    self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full()
        if cpcheck:
            self.check_cp(chemical_potential)
        comm.barrier()

    def test_pm(self):
        '''Test our ability to compute the density matrix with PM.'''
        self.basic_solver(nt.DensityMatrixSolvers.PM)

    def test_trs2(self):
        '''Test routines to compute the density matrix with TRS2.'''
        self.basic_solver(nt.DensityMatrixSolvers.TRS2)

    def test_trs4(self):
        '''Test routines to compute the density matrix with TRS4.'''
        self.basic_solver(nt.DensityMatrixSolvers.TRS4)

    def test_HPCP(self):
        '''Test routines to compute the density matrix with HPCP.'''
        self.basic_solver(nt.DensityMatrixSolvers.HPCP)

    def test_dense(self):
        '''Test routines to compute the density matrix with eigensolver.'''
        self.basic_solver(nt.DensityMatrixSolvers.DenseSolver)

    def test_energy_density(self):
        '''Test the routines to compute the weighted-energy density matrix.'''
        # Reference Solution
        self.create_matrices()
        fmat = mmread(self.hamiltonian)
        dmat = mmread(self.density)
        edm = dmat.dot(fmat).dot(dmat)

        # NTPoly
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        density_matrix = nt.Matrix_ps(self.density)
        edm_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        nt.DensityMatrixSolvers.EnergyDensityMatrix(
            fock_matrix, density_matrix, edm_matrix)

        edm_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        # Compare
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            normval = abs(norm(edm - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        comm.barrier()


class TestChemistry_r(TestChemistry, unittest.TestCase):
    '''Specialization for real matrices'''
    # complex test
    is_complex = False

    def testrealio(self):
        '''Test routines to read data produced by a real chemistry program.'''
        density_matrix = nt.Matrix_ps(self.realio)
        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            self.CheckMat = mmread(self.realio)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        comm.barrier()

    def test_PExtrapolate(self):
        '''Test the density extrapolation routine.'''
        f1 = nt.Matrix_ps(self.geomh1)
        o1 = nt.Matrix_ps(self.geomo1)
        o2 = nt.Matrix_ps(self.geomo2)
        isqm1 = nt.Matrix_ps(f1.GetActualDimension())
        d1 = nt.Matrix_ps(f1.GetActualDimension())
        extrapd = nt.Matrix_ps(f1.GetActualDimension())

        permutation = nt.Permutation(f1.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(
            o1, isqm1, self.solver_parameters)
        nt.DensityMatrixSolvers.TRS2(
            f1, isqm1, self.nel, d1, self.solver_parameters)

        nt.GeometryOptimization.PurificationExtrapolate(d1, o2, self.nel,
                                                        extrapd,
                                                        self.solver_parameters)
        extrapd.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full_extrap()
        comm.barrier()

    def test_SExtrapolate(self):
        '''Test the density extrapolation routine.'''
        f1 = nt.Matrix_ps(self.geomh1)
        o1 = nt.Matrix_ps(self.geomo1)
        o2 = nt.Matrix_ps(self.geomo2)
        isqm1 = nt.Matrix_ps(f1.GetActualDimension())
        d1 = nt.Matrix_ps(f1.GetActualDimension())
        extrapd = nt.Matrix_ps(f1.GetActualDimension())

        permutation = nt.Permutation(f1.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(
            o1, isqm1, self.solver_parameters)
        nt.DensityMatrixSolvers.TRS2(
            f1, isqm1, self.nel, d1, self.solver_parameters)

        nt.GeometryOptimization.LowdinExtrapolate(d1, o1, o2, extrapd,
                                                  self.solver_parameters)
        extrapd.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full_extrap()
        comm.barrier()


class TestChemistry_c(TestChemistry, unittest.TestCase):
    '''Specialization for complex matrices.'''
    is_complex = True


if __name__ == '__main__':
    unittest.main()
