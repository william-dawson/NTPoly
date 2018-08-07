'''@package testDistributedSparseMatrix
A test suite for the Distributed Sparse Matrix module.'''
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
# MPI global communicator.
comm = MPI.COMM_WORLD

from Helpers import THRESHOLD
from Helpers import result_file
from Helpers import scratch_dir


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


class TestDistributedMatrix(unittest.TestCase):
    '''A test class for the distributed matrix module.'''
    # Parameters for the tests
    parameters = []
    # Where to store the result file
    result_file = scratch_dir + "/result.mtx"
    # Matrix to compare against
    CheckMat = 0
    # Rank of the current process.
    my_rank = 0

    @classmethod
    def setUpClass(self):
        '''Set up tests.'''
        self.process_rows = int(os.environ['PROCESS_ROWS'])
        self.process_columns = int(os.environ['PROCESS_COLUMNS'])
        self.process_slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(
            self.process_rows, self.process_columns, self.process_slices)
        # Make sure we can destruct without any problems.
        nt.DestructProcessGrid()
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
            if (self.my_rank == 0):
                matrix1 = random(param.rows, param.columns,
                                 param.sparsity, format="csr")
                mmwrite(scratch_dir + "/matrix1.mtx", csr_matrix(matrix1))
                self.CheckMat = matrix1

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix1.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_readwritebinary(self):
        '''Test our ability to read and write binary.'''
        for param in self.parameters:
            if (self.my_rank == 0):
                matrix1 = random(param.rows, param.columns,
                                 param.sparsity, format="csr")
                mmwrite(scratch_dir + "/matrix1.mtx",
                        csr_matrix(matrix1), symmetry="general")
                self.CheckMat = matrix1

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix1.WriteToBinary(scratch_dir + "/matrix2.mtx")
            ntmatrix2 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix2.mtx", True)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    def test_gettripletlist(self):
        '''Test extraction of triplet list.'''
        for param in self.parameters:
            if (self.my_rank == 0):
                matrix1 = random(param.rows, param.columns,
                                 param.sparsity, format="csr")
                mmwrite(scratch_dir + "/matrix1.mtx", csr_matrix(matrix1))
                self.CheckMat = matrix1

            comm.barrier()
            if param.sparsity > 0.0:
                ntmatrix1 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix1.mtx", False)
            else:
                ntmatrix1 = nt.DistributedSparseMatrix(param.rows)
            triplet_list = nt.TripletList(0)
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
            if (self.my_rank == 0):
                matrix1 = random(param.rows, param.columns,
                                 param.sparsity, format="csr")
                mmwrite(scratch_dir + "/matrix1.mtx", csr_matrix(matrix1))
                self.CheckMat = matrix1

            comm.barrier()

            if param.sparsity > 0.0:
                ntmatrix1 = nt.DistributedSparseMatrix(
                    scratch_dir + "/matrix1.mtx", False)
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

            triplet_list = nt.TripletList(0)
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
            if (self.my_rank == 0):
                matrix1 = random(param.rows, param.columns,
                                 param.sparsity, format="csr")
                mmwrite(scratch_dir + "/matrix1.mtx", csr_matrix(matrix1))
                self.CheckMat = matrix1.T

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix2.Transpose(ntmatrix1)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()


if __name__ == '__main__':
    unittest.main()
