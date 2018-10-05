'''
@package test_matrix
A test suite for paralle matrices.
'''
import unittest
import NTPolySwig as nt
from random import randrange, seed, sample
import scipy
import scipy.sparse
from scipy.sparse import random, csr_matrix
from scipy.sparse.linalg import norm
from scipy.io import mmread, mmwrite
from numpy import zeros
import os
import sys
from mpi4py import MPI
from helpers import THRESHOLD
from helpers import result_file
from helpers import scratch_dir
# MPI global communicator.
comm = MPI.COMM_WORLD


class TestParameters:
    '''An internal class for holding internal class parameters.'''

    def __init__(self, rows, columns, sparsity):
        '''Default constructor.
        @param[in] rows matrix rows.
        @param[in] columns matrix columns.
        @param[in] sparsity matrix sparsity.
        '''
        self.rows = rows
        self.columns = columns
        self.sparsity = sparsity

    def create_matrix(self, complex=False):
        '''
        Create the test matrix with the following parameters.
        '''
        r = self.rows
        c = self.columns
        s = self.sparsity
        if complex:
            mat = random(r, c, s, format="csr")
            mat += 1j * random(r, c, s, format="csr")
        else:
            mat = random(r, c, s, format="csr")
        return csr_matrix(mat)


class TestPSMatrix(unittest.TestCase):
    '''A test class for parallel matrices.'''
    # Parameters for the tests
    parameters = []
    # Input file name 1
    input_file1 = scratch_dir + "/matrix1.mtx"
    # Input file name 2
    input_file2 = scratch_dir + "/matrix2.mtx"
    # Input file name 3
    input_file3 = scratch_dir + "/matrix3.mtx"
    # Where to store the result file
    result_file = scratch_dir + "/result.mtx"
    # Matrix to compare against
    CheckMat = 0
    # Rank of the current process.
    my_rank = 0
    # Type of triplets to use
    TripletList = nt.TripletList_r
    # Whether the matrix is complex or not
    complex = False

    def write_matrix(self, mat, file_name):
        if self.my_rank == 0:
            mmwrite(file_name, csr_matrix(mat))
        comm.barrier()

    @classmethod
    def setUpClass(self):
        '''Set up test suite.'''
        rows = int(os.environ['PROCESS_ROWS'])
        columns = int(os.environ['PROCESS_COLUMNS'])
        slices = int(os.environ['PROCESS_SLICES'])
        # global process grid
        nt.ConstructGlobalProcessGrid(rows, columns, slices)

    @classmethod
    def tearDownClass(self):
        '''Cleanup this test'''
        nt.DestructGlobalProcessGrid()

    def setUp(self):
        '''Set up specific tests.'''
        mat_size = 33
        self.process_rows = int(os.environ['PROCESS_ROWS'])
        self.process_columns = int(os.environ['PROCESS_COLUMNS'])
        self.process_slices = int(os.environ['PROCESS_SLICES'])

        self.grid = nt.ProcessGrid(
            self.process_rows, self.process_columns, self.process_slices)
        self.myrow = self.grid.GetMyRow()
        self.mycolumn = self.grid.GetMyColumn()
        self.myslice = self.grid.GetMySlice()

        self.my_rank = comm.Get_rank()
        self.parameters = []
        self.parameters.append(TestParameters(mat_size, mat_size, 1.0))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.2))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.0))

    def tearDown(self):
        '''Cleanup this test.'''
        del self.grid

    def check_result(self):
        '''Compare two matrices.'''
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(self.result_file)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    def test_grid(self):
        '''Test the simplified process grid interface'''
        self.assertEqual(self.process_rows, nt.GetGlobalNumRows())
        self.assertEqual(self.process_columns, nt.GetGlobalNumColumns())
        self.assertEqual(self.process_slices, nt.GetGlobalNumSlices())

        self.assertEqual(self.process_rows, self.grid.GetNumRows())
        self.assertEqual(self.process_columns, self.grid.GetNumColumns())
        self.assertEqual(self.process_slices, self.grid.GetNumSlices())

        total_procs = self.process_rows * self.process_columns * \
            self.process_slices
        new_grid = nt.ProcessGrid(self.process_slices)
        new_total_procs = new_grid.GetNumRows() * new_grid.GetNumColumns() * \
            new_grid.GetNumSlices()
        self.assertEqual(total_procs, new_total_procs)
        del new_grid

    def test_read(self):
        '''Test our ability to read and write matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            ntmatrix1.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_read_pg(self):
        '''Test our ability to read and write matrices on a given grid.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            ntmatrix1 = nt.Matrix_ps(self.input_file1, self.grid, False)
            ntmatrix1.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_copy_grid(self):
        '''Test process grid copying'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            new_grid = nt.ProcessGrid(self.grid)

            ntmatrix1 = nt.Matrix_ps(self.input_file1, new_grid, False)
            ntmatrix1.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_readwritebinary(self):
        '''Test our ability to read and write binary.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            ntmatrix1.WriteToBinary(self.input_file2)
            ntmatrix2 = nt.Matrix_ps(self.input_file2, True)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_readwritebinary_pg(self):
        '''Test our ability to read and write binary on a given grid.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            ntmatrix1 = nt.Matrix_ps(self.input_file1, self.grid, False)
            ntmatrix1.WriteToBinary(self.input_file2)
            ntmatrix2 = nt.Matrix_ps(self.input_file2, self.grid, True)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_gettripletlist(self):
        '''Test extraction of triplet list.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)

            triplet_list = self.TripletList(0)
            if self.myslice == 0:
                ntmatrix1.GetTripletList(triplet_list)
            ntmatrix2 = nt.Matrix_ps(ntmatrix1.GetActualDimension())
            ntmatrix2.FillFromTripletList(triplet_list)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_repartition(self):
        '''Test extraction of triplet list via repartition function.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)

            # Compute a random permutation
            seed_val = randrange(sys.maxsize)
            global_seed = comm.bcast(seed_val, root=0)
            seed(global_seed)
            dimension = ntmatrix1.GetActualDimension()
            row_end_list = sample(range(0, dimension - 1),
                                  self.process_rows - 1)
            col_end_list = sample(range(0, dimension - 1),
                                  self.process_columns - 1)
            row_end_list.append(dimension)
            col_end_list.append(dimension)
            row_start_list = [0]
            for i in range(1, len(row_end_list)):
                row_start_list.append(row_end_list[i - 1])
            col_start_list = [0]
            for i in range(1, len(col_end_list)):
                col_start_list.append(col_end_list[i - 1])

            triplet_list = self.TripletList(0)
            if self.myslice == 0:
                ntmatrix1.GetMatrixBlock(triplet_list,
                                         row_start_list[self.myrow],
                                         row_end_list[self.myrow],
                                         col_start_list[self.mycolumn],
                                         col_end_list[self.mycolumn])

            ntmatrix2 = nt.Matrix_ps(ntmatrix1.GetActualDimension())
            ntmatrix2.FillFromTripletList(triplet_list)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_slice(self):
        '''Test slicing of a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            if param.sparsity > 0.0:
                ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            else:
                ntmatrix1 = nt.Matrix_ps(param.rows)

            # Compute a random slicing
            seed_val = randrange(sys.maxsize)
            global_seed = comm.bcast(seed_val, root=0)
            seed(global_seed)
            dimension = ntmatrix1.GetActualDimension()
            end_row = sample(range(1, dimension - 1), 1)[0]
            end_col = sample(range(1, dimension - 1), 1)[0]
            start_row = sample(range(0, end_row), 1)[0]
            start_col = sample(range(0, end_col), 1)[0]

            # Compute the reference result
            sub_mat = matrix1[start_row:end_row + 1, start_col:end_col + 1]
            new_dim = max(end_row - start_row + 1, end_col - start_col + 1)
            space_mat = zeros((new_dim, new_dim))
            if self.complex:
                space_mat = 1j * space_mat
            space_mat[:end_row - start_row + 1, :end_col -
                      start_col + 1] = sub_mat.todense()
            self.CheckMat = csr_matrix(space_mat)

            # Compute with ntpoly
            ntmatrix2 = nt.Matrix_ps(ntmatrix1.GetActualDimension())
            ntmatrix1.GetMatrixSlice(
                ntmatrix2, start_row, end_row, start_col, end_col)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_transpose(self):
        '''Test our ability to transpose matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)

            self.CheckMat = matrix1.T
            ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            ntmatrix2 = nt.Matrix_ps(ntmatrix1.GetActualDimension())
            ntmatrix2.Transpose(ntmatrix1)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()


class TestPSMatrix_c(TestPSMatrix):
    '''Specialization for complex matrices'''
    TripletList = nt.TripletList_c
    complex = True

    def test_conjugatetranspose(self):
        '''Test our ability to compute the conjugate transpose of a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)

            self.CheckMat = matrix1.H

            ntmatrix1 = nt.Matrix_ps(self.input_file1, False)
            ntmatrix2 = nt.Matrix_ps(ntmatrix1.GetActualDimension())
            ntmatrix2.Transpose(ntmatrix1)
            ntmatrix2.Conjugate()
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()


if __name__ == '__main__':
    unittest.main()
