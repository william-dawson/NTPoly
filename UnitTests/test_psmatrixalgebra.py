'''
@package test_psmatrixalgebra
A test suite for parallel matrix algebra.
'''
from helpers import THRESHOLD
from helpers import result_file
from helpers import scratch_dir
import unittest
import NTPolySwig as nt
import scipy
from scipy.sparse import random, csr_matrix
from scipy.sparse.linalg import norm
from scipy.io import mmread, mmwrite
from numpy import sum, multiply, conj
import os
from mpi4py import MPI
# MPI global communicator.
comm = MPI.COMM_WORLD


rows = int(os.environ['PROCESS_ROWS'])
columns = int(os.environ['PROCESS_COLUMNS'])
slices = int(os.environ['PROCESS_SLICES'])

class TestParameters:
    '''An internal class for holding internal class parameters.'''

    def __init__(self, rows, sparsity, sparsity2):
        '''Default constructor.
        @param[in] rows matrix rows.
        @param[in] columns matrix columns.
        @param[in] sparsity matrix sparsity.
        '''
        self.rows = rows
        self.columns = rows
        self.sparsity = sparsity
        self.sparsity2 = sparsity2

    def create_matrix(self, snum=1, complex=False):
        '''
        Create the test matrix with the following parameters.
        '''
        r = self.rows
        c = self.columns
        if snum == 1:
            s = self.sparsity
        else:
            s = self.sparsity2
        if complex:
            mat = random(r, c, s, format="csr")
            mat += 1j * random(r, c, s, format="csr")
        else:
            mat = random(r, c, s, format="csr")

        return mat


class TestPSMatrixAlgebra:
    '''A test class for parallel matrix algebra.'''
    # Parameters for the tests
    parameters = []
    # Place to store the result matrix.
    result_file = scratch_dir + "/result.mtx"
    # Input file name 1
    input_file1 = scratch_dir + "/matrix1.mtx"
    # Input file name 2
    input_file2 = scratch_dir + "/matrix2.mtx"
    # Input file name 3
    input_file3 = scratch_dir + "/matrix3.mtx"
    # Matrix to compare against.
    CheckMat = 0
    # Rank of the current process.
    my_rank = 0
    # Whether the first matrix is complex or not
    complex1 = False
    # Whether the second matrix is complex or not
    complex2 = False
    # The size of the matrix
    mat_size = 33

    def write_matrix(self, mat, file_name):
        if self.my_rank == 0:
            mmwrite(file_name, csr_matrix(mat))
        comm.barrier()

    @classmethod
    def setUpClass(self):
        '''Set up test suite.'''
        nt.ConstructGlobalProcessGrid(rows, columns, slices)

    @classmethod
    def tearDownClass(self):
        '''Cleanup this test suite.'''
        nt.DestructGlobalProcessGrid()

    def tearDown(self):
        '''Cleanup this test.'''
        del self.grid

    def setUp(self):
        '''Set up a specific test.'''
        self.grid = nt.ProcessGrid(rows, columns, slices)
        self.my_rank = comm.Get_rank()

        self.parameters = []
        self.parameters.append(TestParameters(self.mat_size, 1.0, 1.0))
        self.parameters.append(TestParameters(self.mat_size, 0.2, 0.2))
        self.parameters.append(TestParameters(self.mat_size, 0.0, 0.0))
        self.parameters.append(TestParameters(self.mat_size, 1.0, 0.0))
        self.parameters.append(TestParameters(self.mat_size, 0.0, 1.0))

    def check_result(self):
        '''Compare two matrices.'''
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(self.result_file)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    def test_addition_pg(self):
        '''Test routines to add together matrices with an explicit grid.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(snum=1, complex=self.complex1)
            matrix2 = param.create_matrix(snum=2, complex=self.complex2)
            self.write_matrix(matrix1, self.input_file1)
            self.write_matrix(matrix2, self.input_file2)

            self.CheckMat = matrix1 + matrix2

            comm.barrier()
            ntmatrix1 = nt.Matrix_ps(self.input_file1, self.grid, False)
            ntmatrix2 = nt.Matrix_ps(self.input_file2, self.grid, False)
            ntmatrix2.Increment(ntmatrix1)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_addition(self):
        '''Test routines to add together matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(snum=1, complex=self.complex1)
            matrix2 = param.create_matrix(snum=2, complex=self.complex2)
            self.write_matrix(matrix1, self.input_file1)
            self.write_matrix(matrix2, self.input_file2)

            self.CheckMat = matrix1 + matrix2

            comm.barrier()
            ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            ntmatrix2 = nt.Matrix_ps(self.input_file2, False)
            ntmatrix2.Increment(ntmatrix1)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_pairwisemultiply(self):
        '''Test routines to pairwise multiply two matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(snum=1, complex=self.complex1)
            matrix2 = param.create_matrix(snum=2, complex=self.complex2)
            self.write_matrix(matrix1, self.input_file1)
            self.write_matrix(matrix2, self.input_file2)
            self.CheckMat = csr_matrix(multiply(
                matrix1.todense(), matrix2.todense()))

            comm.barrier()
            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)
            if param.sparsity2 > 0.0:
                ntmatrix2 = nt.Matrix_ps(self.input_file2, False)
            else:
                ntmatrix2 = nt.Matrix_ps(param.rows)
            ntmatrix3 = nt.Matrix_ps(param.rows)
            ntmatrix3.PairwiseMultiply(ntmatrix1, ntmatrix2)
            ntmatrix3.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_multiply(self):
        '''Test routines to multiply two matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(snum=1, complex=self.complex1)
            matrix2 = param.create_matrix(snum=2, complex=self.complex2)
            self.write_matrix(matrix1, self.input_file1)
            self.write_matrix(matrix2, self.input_file2)

            self.CheckMat = matrix1.dot(matrix2)
            comm.barrier()

            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)
            if param.sparsity2 > 0.0:
                ntmatrix2 = nt.Matrix_ps(self.input_file2, False)
            else:
                ntmatrix2 = nt.Matrix_ps(param.rows)
            ntmatrix3 = nt.Matrix_ps(param.rows)
            memory_pool = nt.PMatrixMemoryPool(ntmatrix1)
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_reverse(self):
        '''Test routines to permute a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(snum=1, complex=self.complex1)
            self.write_matrix(matrix1, self.input_file1)

            self.CheckMat = matrix1
            comm.barrier()

            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)
            permute_rows = nt.Matrix_ps(param.rows)
            permute_columns = nt.Matrix_ps(param.rows)
            temp_matrix = nt.Matrix_ps(param.rows)
            memory_pool = nt.PMatrixMemoryPool(ntmatrix1)
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


class TestPSMatrixAlgebra_r(TestPSMatrixAlgebra, unittest.TestCase):
    '''Special routines for real algebra'''
    # Whether the first matrix is complex or not
    complex1 = False
    # Whether the second matrix is complex or not
    complex2 = False

    def test_dot(self):
        '''Test routines to add together matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(snum=1, complex=self.complex1)
            matrix2 = param.create_matrix(snum=2, complex=self.complex2)
            self.write_matrix(matrix1, self.input_file1)
            self.write_matrix(matrix2, self.input_file2)
            check = 0
            if self.my_rank == 0:
                check = sum(multiply(matrix1.todense(), matrix2.todense()))

            check = comm.bcast(check, root=0)
            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)
            if param.sparsity2 > 0.0:
                ntmatrix2 = nt.Matrix_ps(self.input_file2, False)
            else:
                ntmatrix2 = nt.Matrix_ps(param.rows)

            result = 0
            result = ntmatrix2.Dot(ntmatrix1)
            comm.barrier()

            normval = abs(result - check)

            self.assertLessEqual(normval, THRESHOLD)


class TestPSMatrixAlgebra_c(TestPSMatrixAlgebra, unittest.TestCase):
    '''Specialization for complex algebra'''
    # Whether the first matrix is complex or not
    complex1 = True
    # Whether the second matrix is complex or not
    complex2 = True

    def test_dot(self):
        '''Test routines to add together matrices.'''
        for param in self.parameters:
            comm.barrier()
            matrix1 = param.create_matrix(snum=1, complex=self.complex1)
            matrix2 = param.create_matrix(snum=2, complex=self.complex2)
            self.write_matrix(matrix1, self.input_file1)
            self.write_matrix(matrix2, self.input_file2)

            check = 0
            if self.my_rank == 0:
                check = sum(
                    multiply(conj(matrix1.todense()), matrix2.todense()))

            comm.barrier()
            check = comm.bcast(check, root=0)
            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)
            if param.sparsity2 > 0.0:
                ntmatrix2 = nt.Matrix_ps(self.input_file2, False)
            else:
                ntmatrix2 = nt.Matrix_ps(param.rows)

            result = ntmatrix2.Dot(ntmatrix1)
            comm.barrier()

            check = check.real

            normval = abs(result - check)
            self.assertLessEqual(normval, THRESHOLD)


class TestPSMatrixAlgebra_rc(TestPSMatrixAlgebra, unittest.TestCase):
    '''Specialization for real-complex mixing'''
    # Whether the first matrix is complex or not
    complex1 = True
    # Whether the second matrix is complex or not
    complex2 = False

    def setUp(self):
        '''Set up a specific test.'''
        self.grid = nt.ProcessGrid(rows, columns, slices)
        self.my_rank = comm.Get_rank()

        self.parameters = [TestParameters(self.mat_size, 0.2, 0.2)]


class TestPSMatrixAlgebra_cr(TestPSMatrixAlgebra, unittest.TestCase):
    '''Specialization for complex-real mixing'''
    # Whether the first matrix is complex or not
    complex1 = False
    # Whether the second matrix is complex or not
    complex2 = True

    def setUp(self):
        '''Set up a specific test.'''
        self.grid = nt.ProcessGrid(rows, columns, slices)
        self.my_rank = comm.Get_rank()

        self.parameters = []
        self.parameters = [TestParameters(self.mat_size, 0.2, 0.2)]


if __name__ == '__main__':
    unittest.main()
