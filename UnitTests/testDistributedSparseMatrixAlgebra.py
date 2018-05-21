'''@package testDistributedSparseMatrix
A test suite for the Distributed Sparse Matrix module.'''
from Helpers import THRESHOLD
from Helpers import result_file
from Helpers import scratch_dir
import unittest
import NTPolySwig as nt
import scipy
from scipy.sparse import random, csr_matrix
from scipy.sparse.linalg import norm
from scipy.io import mmread, mmwrite
from numpy import sum, multiply
import os
from mpi4py import MPI
# MPI global communicator.
comm = MPI.COMM_WORLD


class TestParameters:
    '''An internal class for holding internal class parameters.'''

    def __init__(self, rows, columns, sparsity, sparsity2):
        '''Default constructor.
        @param[in] rows matrix rows.
        @param[in] columns matrix columns.
        @param[in] sparsity matrix sparsity.
        '''
        self.rows = rows
        self.columns = columns
        self.sparsity = sparsity
        self.sparsity2 = sparsity2


class TestDistributedMatrixAlgebra(unittest.TestCase):
    '''A test class for the distributed matrix module.'''
    # Parameters for the tests
    parameters = []
    # Place to store the result matrix.
    result_file = scratch_dir + "/result.mtx"
    # Matrix to compare against.
    CheckMat = 0
    # Rank of the current process.
    my_rank = 0

    @classmethod
    def setUpClass(self):
        '''Set up tests.'''
        rows = int(os.environ['PROCESS_ROWS'])
        columns = int(os.environ['PROCESS_COLUMNS'])
        slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(rows, columns, slices)

    def setUp(self):
        '''Set up a specific test.'''
        mat_size = 64
        self.my_rank = comm.Get_rank()
        self.parameters.append(TestParameters(mat_size, mat_size, 1.0, 1.0))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.2, 0.2))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.0, 0.0))
        self.parameters.append(TestParameters(mat_size, mat_size, 1.0, 0.0))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.0, 1.0))
        self.parameters.append(TestParameters(7, 7, 0.2, 0.2))

    def check_result(self):
        '''Compare two matrices.'''
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(self.result_file)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    def test_addition(self):
        '''Test our ability to add together matrices.'''
        for param in self.parameters:
            matrix1 = random(param.rows, param.columns,
                             param.sparsity, format="csr")
            matrix2 = random(param.rows, param.columns,
                             param.sparsity2, format="csr")
            self.CheckMat = matrix1 + matrix2
            if self.my_rank == 0:
                mmwrite(scratch_dir + "/matrix1.mtx",
                        csr_matrix(matrix1), symmetry="general")
                mmwrite(scratch_dir + "/matrix2.mtx",
                        csr_matrix(matrix2), symmetry="general")

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix2.mtx", False)
            ntmatrix2.Increment(ntmatrix1)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_dot(self):
        '''Test our ability to add together matrices.'''
        for param in self.parameters:
            comm.barrier()
            check = 0
            if self.my_rank == 0:
                matrix1 = random(param.rows, param.columns,
                                 param.sparsity, format="csr")
                matrix2 = random(param.rows, param.columns,
                                 param.sparsity2, format="csr")
                check = sum(multiply(matrix1.todense(), matrix2.todense()))
                mmwrite(scratch_dir + "/matrix1.mtx",
                        csr_matrix(matrix1), symmetry="general")
                mmwrite(scratch_dir + "/matrix2.mtx",
                        csr_matrix(matrix2), symmetry="general")

            comm.barrier()
            check = comm.bcast(check, root=0)
            if param.sparsity > 0.0:
                ntmatrix1 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix1.mtx", False)
            else:
                ntmatrix1 = nt.DistributedSparseMatrix(param.rows)
            if param.sparsity2 > 0.0:
                ntmatrix2 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix2.mtx", False)
            else:
                ntmatrix2 = nt.DistributedSparseMatrix(param.rows)

            result = ntmatrix2.Dot(ntmatrix1)
            comm.barrier()

            normval = abs(result - check)

            self.assertLessEqual(normval, THRESHOLD)

    def test_pairwisemultiply(self):
        '''Test our ability to pairwise multiply two matrices.'''
        for param in self.parameters:
            matrix1 = random(param.rows, param.columns,
                             param.sparsity, format="csr")
            matrix2 = random(param.columns, param.rows,
                             param.sparsity2, format="csr")
            self.CheckMat = csr_matrix(multiply(
                matrix1.todense(), matrix2.todense()))
            if self.my_rank == 0:
                mmwrite(scratch_dir + "/matrix1.mtx",
                        csr_matrix(matrix1), symmetry="general")
                mmwrite(scratch_dir + "/matrix2.mtx",
                        csr_matrix(matrix2), symmetry="general")

            comm.barrier()
            if param.sparsity > 0.0:
                ntmatrix1 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix1.mtx", False)
            else:
                ntmatrix1 = nt.DistributedSparseMatrix(param.rows)
            if param.sparsity2 > 0.0:
                ntmatrix2 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix2.mtx", False)
            else:
                ntmatrix2 = nt.DistributedSparseMatrix(param.rows)
            ntmatrix3 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix3.PairwiseMultiply(ntmatrix1, ntmatrix2)
            ntmatrix3.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_multiply(self):
        '''Test our ability to multiply two matrices.'''
        for param in self.parameters:
            matrix1 = random(param.rows, param.columns,
                             param.sparsity, format="csr")
            matrix2 = random(param.columns, param.rows,
                             param.sparsity2, format="csr")
            self.CheckMat = matrix1.dot(matrix2)
            if self.my_rank == 0:
                mmwrite(scratch_dir + "/matrix1.mtx",
                        csr_matrix(matrix1), symmetry="general")
                mmwrite(scratch_dir + "/matrix2.mtx",
                        csr_matrix(matrix2), symmetry="general")

            comm.barrier()
            if param.sparsity > 0.0:
                ntmatrix1 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix1.mtx", False)
            else:
                ntmatrix1 = nt.DistributedSparseMatrix(param.rows)
            if param.sparsity2 > 0.0:
                ntmatrix2 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix2.mtx", False)
            else:
                ntmatrix2 = nt.DistributedSparseMatrix(param.rows)
            ntmatrix3 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            memory_pool = nt.DistributedMatrixMemoryPool(ntmatrix1)
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_reverse(self):
        '''Test our ability to permute a matrix.'''
        for param in self.parameters:
            matrix1 = random(param.rows, param.columns,
                             param.sparsity, format="csr")
            self.CheckMat = matrix1

            if self.my_rank == 0:
                mmwrite(scratch_dir + "/matrix1.mtx",
                        csr_matrix(matrix1), symmetry="general")

            comm.barrier()
            if param.sparsity > 0.0:
                ntmatrix1 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix1.mtx", False)
            else:
                ntmatrix1 = nt.DistributedSparseMatrix(param.rows)
            permute_rows = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            permute_columns = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            temp_matrix = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            memory_pool = nt.DistributedMatrixMemoryPool(ntmatrix1)
            permutation = nt.Permutation(ntmatrix1.GetLogicalDimension())

            permutation.SetReversePermutation()
            permute_rows.FillDistributedPermutation(
                permutation, permuterows=True)
            permute_columns.FillDistributedPermutation(
                permutation, permuterows=False)

            temp_matrix.Gemm(permute_rows, ntmatrix1, memory_pool)
            ntmatrix1.Gemm(temp_matrix, permute_columns, memory_pool)
            temp_matrix.Gemm(permute_columns, ntmatrix1, memory_pool)
            ntmatrix1.Gemm(temp_matrix, permute_rows, memory_pool)

            ntmatrix1.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()


if __name__ == '__main__':
    unittest.main()
