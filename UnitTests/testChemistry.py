''' @package testDistributedSparseMatrix
A test suite for the Distributed Sparse Matrix module.
'''
import unittest
import NTPolySwig as nt

import scipy
from scipy.sparse.linalg import norm
from scipy.io import mmread
from scipy.linalg import eigh
import os
from mpi4py import MPI
# MPI globa communicator
comm = MPI.COMM_WORLD

from Helpers import THRESHOLD
from Helpers import result_file
from Helpers import scratch_dir


class TestChemistry(unittest.TestCase):
    '''A test class for the distributed matrix module.'''
    # Parameters for the tests
    parameters = []
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
        self.solver_parameters = nt.IterativeSolverParameters()
        self.solver_parameters.SetVerbosity(True)
        self.hamiltonian = os.environ["HAMILTONIAN"]
        self.overlap = os.environ["OVERLAP"]
        self.density = os.environ["DENSITY"]
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

    def test_trs2(self):
        '''Test our ability to compute the density matrix with TRS2.'''
        fock_matrix = nt.DistributedSparseMatrix(self.hamiltonian)
        overlap_matrix = nt.DistributedSparseMatrix(self.overlap)
        inverse_sqrt_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())
        density_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        chemical_potential = nt.DensityMatrixSolvers.TRS2(fock_matrix,
                                                          inverse_sqrt_matrix,
                                                          self.nel, density_matrix,
                                                          self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full()
        self.check_cp(chemical_potential)

    def test_trs4(self):
        '''Test our ability to compute the density matrix with TRS4.'''
        fock_matrix = nt.DistributedSparseMatrix(self.hamiltonian)
        overlap_matrix = nt.DistributedSparseMatrix(self.overlap)
        inverse_sqrt_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())
        density_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        chemical_potential = nt.DensityMatrixSolvers.TRS4(fock_matrix,
                                                          inverse_sqrt_matrix,
                                                          self.nel, density_matrix,
                                                          self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full()
        self.check_cp(chemical_potential)

    def test_HPCP(self):
        '''Test our ability to compute the density matrix with HPCP.'''
        fock_matrix = nt.DistributedSparseMatrix(self.hamiltonian)
        overlap_matrix = nt.DistributedSparseMatrix(self.overlap)
        inverse_sqrt_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())
        density_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        chemical_potential = nt.DensityMatrixSolvers.HPCP(fock_matrix,
                                                          inverse_sqrt_matrix,
                                                          self.nel, density_matrix,
                                                          self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full()
        self.check_cp(chemical_potential)

    # def test_HPCPPlus(self):
    #     '''Test our ability to compute the density matrix with HPCP+.'''
    #     fock_matrix = nt.DistributedSparseMatrix(self.hamiltonian)
    #     overlap_matrix = nt.DistributedSparseMatrix(self.overlap)
    #     inverse_sqrt_matrix = nt.DistributedSparseMatrix(
    #         fock_matrix.GetActualDimension())
    #     density_matrix = nt.DistributedSparseMatrix(
    #         fock_matrix.GetActualDimension())
    #
    #     permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
    #     permutation.SetRandomPermutation()
    #     self.solver_parameters.SetLoadBalance(permutation)
    #
    #     nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
    #                                            inverse_sqrt_matrix,
    #                                            self.solver_parameters)
    #     chemical_potential = nt.DensityMatrixSolvers.HPCPPlus(fock_matrix,
    #                                                           inverse_sqrt_matrix,
    #                                                           self.nel, density_matrix,
    #                                                           self.solver_parameters)
    #
    #     density_matrix.WriteToMatrixMarket(result_file)
    #     comm.barrier()
    #
    #     self.check_full()
    #     self.check_cp(chemical_potential)

    def test_cg(self):
        '''Test our ability to compute the density matrix with conjugate
           gradient.'''
        fock_matrix = nt.DistributedSparseMatrix(self.hamiltonian)
        overlap_matrix = nt.DistributedSparseMatrix(self.overlap)
        inverse_sqrt_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())
        density_matrix = nt.DistributedSparseMatrix(
            fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        chemical_potential = nt.MinimizerSolvers.ConjugateGradient(fock_matrix,
                                                                   inverse_sqrt_matrix,
                                                                   self.nel, density_matrix,
                                                                   self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full()
        self.check_cp(chemical_potential)


if __name__ == '__main__':
    unittest.main()
