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
class TestDistributedMatrixAlgebra(unittest.TestCase):
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
        global_norm = comm.bcast(norm, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    ##########################################################################
    # Test our ability to add together matrices.
    #  @param[in] self pointer.

    def test_addition(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            matrix2 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            self.CheckMat = matrix1 + matrix2
            if self.my_rank == 0:
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1),
                                 symmetry="general")
                scipy.io.mmwrite(scratch_dir + "/matrix2.mtx",
                                 scipy.sparse.csr_matrix(matrix2),
                                 symmetry="general")

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix2.mtx", False)
            ntmatrix2.Increment(ntmatrix1)
            ntmatrix2.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()
    ##########################################################################
    # Test our ability to add together matrices.
    #  @param[in] self pointer.

    def test_dot(self):
        for param in self.parameters:
            comm.barrier()
            check = 0
            if self.my_rank == 0:
                matrix1 = scipy.sparse.random(param.rows, param.columns,
                                              param.sparsity,
                                              format="csr")
                matrix2 = scipy.sparse.random(param.rows, param.columns,
                                              param.sparsity,
                                              format="csr")
                check = numpy.sum(numpy.multiply(
                    matrix1.todense(), matrix2.todense()))
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1),
                                 symmetry="general")
                scipy.io.mmwrite(scratch_dir + "/matrix2.mtx",
                                 scipy.sparse.csr_matrix(matrix2),
                                 symmetry="general")

            comm.barrier()
            check = comm.bcast(check, root=0)
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix2.mtx", False)
            result = ntmatrix2.Dot(ntmatrix1)
            comm.barrier()

            norm = abs(result - check)

            self.assertLessEqual(norm, THRESHOLD)
    ##########################################################################
    # Test our ability to pairwise multiply two matrices.
    #  @param[in] self pointer.

    def test_pairwisemultiply(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            matrix2 = scipy.sparse.random(param.columns, param.rows,
                                          param.sparsity,
                                          format="csr")
            self.CheckMat = scipy.sparse.csr_matrix(numpy.multiply(
                matrix1.todense(), matrix2.todense()))
            if self.my_rank == 0:
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1),
                                 symmetry="general")
                scipy.io.mmwrite(scratch_dir + "/matrix2.mtx",
                                 scipy.sparse.csr_matrix(matrix2),
                                 symmetry="general")

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix2.mtx", False)
            ntmatrix3 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            ntmatrix3.PairwiseMultiply(ntmatrix1, ntmatrix2)
            ntmatrix3.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()
    ##########################################################################
    # Test our ability to multiply two matrices.
    #  @param[in] self pointer.

    def test_multiply(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns, param.sparsity,
                                          format="csr")
            matrix2 = scipy.sparse.random(param.columns, param.rows, param.sparsity,
                                          format="csr")
            self.CheckMat = matrix1.dot(matrix2)
            if self.my_rank == 0:
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1), symmetry="general")
                scipy.io.mmwrite(scratch_dir + "/matrix2.mtx",
                                 scipy.sparse.csr_matrix(matrix2), symmetry="general")

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            ntmatrix2 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix2.mtx", False)
            ntmatrix3 = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            memory_pool = nt.DistributedMatrixMemoryPool()
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.result_file)
            comm.barrier()

            self.check_result()

    ##########################################################################
    # Test our ability to permute a matrix.
    #  @param[in] self pointer.
    def test_reverse(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            self.CheckMat = matrix1

            if self.my_rank == 0:
                scipy.io.mmwrite(scratch_dir + "/matrix1.mtx",
                                 scipy.sparse.csr_matrix(matrix1),
                                 symmetry="general")

            comm.barrier()
            ntmatrix1 = nt.DistributedSparseMatrix(
                scratch_dir + "/matrix1.mtx", False)
            permute_rows = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            permute_columns = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            temp_matrix = nt.DistributedSparseMatrix(
                ntmatrix1.GetActualDimension())
            memory_pool = nt.DistributedMatrixMemoryPool()
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
