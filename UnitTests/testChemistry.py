'''
@package testDistributedSparseMatrix
A test suite for the Distributed Sparse Matrix module.
'''
from Helpers import THRESHOLD, EXTRAPTHRESHOLD
from Helpers import result_file
from Helpers import scratch_dir
import unittest
import NTPolySwig as nt
from numpy import diag
import scipy
from scipy.sparse.linalg import norm
from scipy.io import mmread
from scipy.linalg import eigh
import os
from mpi4py import MPI
# MPI global communicator
comm = MPI.COMM_WORLD


class TestChemistry(unittest.TestCase):
    '''A test class for the chemistry solvers.'''
    # Matrix to compare to
    CheckMat = 0
    # Rank of the current process
    my_rank = 0

    @classmethod
    def setUpClass(self):
        '''Set up all of the tests.'''
        rows = int(os.environ['PROCESS_ROWS'])
        columns = int(os.environ['PROCESS_COLUMNS'])
        slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(rows, columns, slices, True)

    def setUp(self):
        '''Set up an individual test.'''
        self.my_rank = comm.Get_rank()
        self.solver_parameters = nt.SolverParameters()
        self.solver_parameters.SetVerbosity(True)
        self.hamiltonian = os.environ["HAMILTONIAN"]
        self.overlap = os.environ["OVERLAP"]
        self.density = os.environ["DENSITY"]
        self.geomh1 = os.environ["GEOMH1"]
        self.geomo1 = os.environ["GEOMO1"]
        self.geomo2 = os.environ["GEOMO2"]
        self.geomd2 = os.environ["GEOMD2"]
        self.realio = os.environ["REALIO"]
        self.nel = 10

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
        pass

    def basic_solver(self, SRoutine):
        '''Test various kinds of density matrix solvers.'''
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
        self.check_cp(chemical_potential)

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

    def test_cg(self):
        '''Test routines to compute the density matrix with conjugate
           gradient.'''
        self.basic_solver(nt.MinimizerSolvers.ConjugateGradient)


class TestChemistry_r(TestChemistry):
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


class TestChemistry_c(TestChemistry):
    '''A test class for complex chemistry solvers.'''

    def setUp(self):
        '''Set up an individual test.'''
        self.my_rank = comm.Get_rank()
        self.solver_parameters = nt.SolverParameters()
        self.solver_parameters.SetVerbosity(True)
        self.hamiltonian = os.environ["HCOMPLEX"]
        self.overlap = os.environ["SCOMPLEX"]
        self.density = os.environ["DCOMPLEX"]
        self.nel = 10


if __name__ == '__main__':
    unittest.main()
