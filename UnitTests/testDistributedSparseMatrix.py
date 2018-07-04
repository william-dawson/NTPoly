'''
@package testDistributedSparseMatrix
A test suite for the Distributed Sparse Matrix module.
'''
import unittest
import NTPolySwig as nt
from random import randrange, seed, sample
import scipy
import scipy.sparse
from scipy.sparse import random, csr_matrix
from scipy.sparse.linalg import norm
from scipy.io import mmread, mmwrite
import os
import sys
from mpi4py import MPI
from Helpers import THRESHOLD
from Helpers import result_file
from Helpers import scratch_dir
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
        r = self.rows
        c = self.columns
        s = self.sparsity
        if complex:
            mat = random(r, c, s, format="csr")
            mat += 1j * random(r, c, s, format="csr")
        else:
            mat = random(r, c, s, format="csr")
        return csr_matrix(mat)


class TestDistributedMatrix(unittest.TestCase):
    '''A test class for the distributed matrix module.'''
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
        '''Set up tests.'''
        self.process_rows = int(os.environ['PROCESS_ROWS'])
        self.process_columns = int(os.environ['PROCESS_COLUMNS'])
        self.process_slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(
            self.process_rows, self.process_columns, self.process_slices)
        self.myrow = nt.GetMyRow()
        self.mycolumn = nt.GetMyColumn()
        self.myslice = nt.GetMySlice()

    def setUp(self):
        '''Set up specific tests.'''
        mat_size = 64
        self.my_rank = comm.Get_rank()
        self.parameters.append(TestParameters(mat_size, mat_size, 1.0))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.2))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.0))
        self.parameters.append(TestParameters(7, 7, 0.2))

    def check_result(self):
        '''Compare two matrices.'''
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(self.result_file)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    def test_read(self):
        '''Test our ability to read and write matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            ntmatrix1 = nt.DistributedSparseMatrix(self.input_file1, False)
            ntmatrix1.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_readwritebinary(self):
        '''Test our ability to read and write binary.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)
            self.CheckMat = matrix1

            ntmatrix1 = nt.DistributedSparseMatrix(self.input_file1, False)
            ntmatrix1.WriteToBinary(self.input_file2)
            ntmatrix2 = nt.DistributedSparseMatrix(self.input_file2, True)
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
                ntmatrix1 = nt.DistributedSparseMatrix(self.input_file1, False)
            else:
                ntmatrix1 = nt.DistributedSparseMatrix(param.rows)

            triplet_list = self.TripletList(0)
            if self.myslice == 0:
                ntmatrix1.GetTripletList(triplet_list)
            ntmatrix2 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
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
                ntmatrix1 = nt.DistributedSparseMatrix(self.input_file1, False)
            else:
                ntmatrix1 = nt.DistributedSparseMatrix(param.rows)

            # Compute a random permutation
            seed_val = randrange(sys.maxsize)
            global_seed = comm.bcast(seed_val, root=0)
            seed(global_seed)
            dimension = ntmatrix1.GetActualDimension()
            row_end_list = sample(range(1, dimension), self.process_rows - 1)
            col_end_list = sample(range(1, dimension),
                                  self.process_columns - 1)
            row_end_list.append(dimension + 1)
            col_end_list.append(dimension + 1)
            row_start_list = [1]
            for i in range(1, len(row_end_list)):
                row_start_list.append(row_end_list[i - 1])
            col_start_list = [1]
            for i in range(1, len(col_end_list)):
                col_start_list.append(col_end_list[i - 1])

            triplet_list = self.TripletList(0)
            if self.myslice == 0:
                ntmatrix1.GetMatrixBlock(triplet_list,
                                         row_start_list[self.myrow],
                                         row_end_list[self.myrow],
                                         col_start_list[self.mycolumn],
                                         col_end_list[self.mycolumn])

            ntmatrix2 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix2.FillFromTripletList(triplet_list)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_transpose(self):
        '''Test our ability to transpose matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)

            self.CheckMat = matrix1.T
            ntmatrix1 = nt.DistributedSparseMatrix(self.input_file1, False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix2.Transpose(ntmatrix1)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()


class TestDistributedMatrix_c(TestDistributedMatrix):
    TripletList = nt.TripletList_c
    complex = True

    def test_conjugatetranspose(self):
        '''Test our ability to compute the conjugate transpose of a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(self.complex)
            self.write_matrix(matrix1, self.input_file1)

            self.CheckMat = matrix1.H

            ntmatrix1 = nt.DistributedSparseMatrix(self.input_file1, False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix2.Transpose(ntmatrix1)
            ntmatrix2.Conjugate()
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()


if __name__ == '__main__':
    unittest.main()
