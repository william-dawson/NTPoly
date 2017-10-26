##########################################################################
# @package testDistributedSparseMatrix
#  A test suite for the Distributed Sparse Matrix module.
import unittest
import NTPolySwig as nt

import scipy
import scipy.sparse
import scipy.io
import time
import numpy
import os
import random
from mpi4py import MPI
comm = MPI.COMM_WORLD

from Helpers import THRESHOLD
from Helpers import result_file
from Helpers import scratch_dir

##########################################################################
# An internal class for holding internal class parameters.


class TestParameters:
    # Default constructor.
    #  @param[in] self pointer.
    #  @param[in] rows matrix rows.
    #  @param[in] columns matrix columns.
    #  @param[in] sparsity matrix sparsity.
    def __init__(self, rows, columns, sparsity):
        self.rows = rows
        self.columns = columns
        self.sparsity = sparsity

##########################################################################
# A test class for the distributed matrix module.


class TestDistributedMatrix(unittest.TestCase):
    # Parameters for the tests
    parameters = []
    result_file = scratch_dir + "/result.mtx"
    CheckMat = 0
    my_rank = 0

    ##########################################################################
    # set up tests
    #  @param[in] self pointer.
    @classmethod
    def setUpClass(self):
        self.process_rows = int(os.environ['PROCESS_ROWS'])
        self.process_columns = int(os.environ['PROCESS_COLUMNS'])
        self.process_slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(
            self.process_rows, self.process_columns, self.process_slices)
        self.myrow = nt.GetMyRow()
        self.mycolumn = nt.GetMyColumn()
        self.myslice = nt.GetMySlice()

    def setUp(self):
        mat_size = 4
        self.my_rank = comm.Get_rank()
        self.parameters.append(TestParameters(mat_size, mat_size, 1.0))
        self.parameters.append(TestParameters(mat_size, mat_size, 0.2))
        self.parameters.append(TestParameters(7, 7, 0.2))

    ##########################################################################
    def check_result(self):
        norm = 0
        if (self.my_rank == 0):
            ResultMat = scipy.io.mmread(self.result_file)
            norm = abs(scipy.sparse.linalg.norm(self.CheckMat - ResultMat))
            #ResultMat = ReadMM(self.result_file)
            #norm = abs(scipy.linalg.norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(norm, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    ##########################################################################
    # Test our ability to read and write matrices.
    #  @param[in] self pointer.
    def test_read(self):
        for param in self.parameters:
            if (self.my_rank == 0):
                matrix1 = scipy.sparse.random(param.rows, param.columns,
                                              param.sparsity,
                                              format="csr")
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1))
                self.CheckMat = matrix1

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix1.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    ##########################################################################
    # Test our ability to read and write binary.
    #  @param[in] self pointer.
    def test_readwritebinary(self):
        for param in self.parameters:
            if (self.my_rank == 0):
                matrix1 = scipy.sparse.random(param.rows, param.columns,
                                              param.sparsity,
                                              format="csr")
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1),
                                 symmetry="general")
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

    ##########################################################################
    # Test extraction of triplet list
    #  @param[in] self pointer.
    def test_gettripletlist(self):
        for param in self.parameters:
            if (self.my_rank == 0):
                matrix1 = scipy.sparse.random(param.rows, param.columns,
                                              param.sparsity,
                                              format="csr")
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1))
                self.CheckMat = matrix1

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            triplet_list = nt.TripletList(0)
            ntmatrix1.GetTripletList(triplet_list)
            ntmatrix2 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix2.FillFromTripletList(triplet_list)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    ##########################################################################
    # Test extraction of triplet list via repartition function
    #  @param[in] self pointer.
    def test_repartition(self):
        for param in self.parameters:
            if (self.my_rank == 0):
                matrix1 = scipy.sparse.random(param.rows, param.columns,
                                              param.sparsity,
                                              format="csr")
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1))
                self.CheckMat = matrix1

            comm.barrier()

            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)

            # Compute a random permutation
            dimension = ntmatrix1.GetActualDimension()
            row_end_list = random.sample(
                range(1, dimension), self.process_rows - 1)
            col_end_list = random.sample(
                range(1, dimension), self.process_columns - 1)
            row_end_list.append(dimension+1)
            col_end_list.append(dimension+1)
            row_start_list = [1]
            for i in range(1, len(row_end_list)):
                row_start_list.append(row_end_list[i-1])
            col_start_list = [1]
            for i in range(1, len(col_end_list)):
                col_start_list.append(col_end_list[i-1])

            triplet_list = nt.TripletList(0)
            ntmatrix1.RepartitionMatrix(triplet_list,
              row_start_list[self.myrow], row_end_list[self.myrow],
              col_start_list[self.mycolumn], col_end_list[self.mycolumn])

            ntmatrix2 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix2.FillFromTripletList(triplet_list)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()
        pass


if __name__ == '__main__':
    unittest.main()
