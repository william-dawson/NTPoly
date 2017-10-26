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
        rows = int(os.environ['PROCESS_ROWS'])
        columns = int(os.environ['PROCESS_COLUMNS'])
        slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(rows, columns, slices)

    def setUp(self):
        mat_size = 64
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
        pass

    ##########################################################################
    # Test extraction of triplet list via repartition function
    #  @param[in] self pointer.
    def test_repartition(self):
        pass
if __name__ == '__main__':
    unittest.main()
