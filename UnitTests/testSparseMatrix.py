##########################################################################
# @package testSparseMatrix
#  A test suite for the Sparse Matrix module.
import unittest
import NTPolySwig as nt

import scipy
import numpy
import scipy.sparse
import scipy.io
import random
import os

from Helpers import THRESHOLD
from Helpers import result_file

##########################################################################
# An internal class for holding test parameters.
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
# A test class for the local matrix module.
class TestLocalMatrix(unittest.TestCase):
    # Parameters for the matrices
    parameters = []
    scratch_dir = os.environ['SCRATCHDIR']
    ##########################################################################
    # set up tests
    #  @param[in] self pointer.

    def setUp(self):
        self.parameters.append(TestParameters(2, 4, 0.0))
        self.parameters.append(TestParameters(8, 8, 0.0))
        self.parameters.append(TestParameters(2, 2, 1.0))
        self.parameters.append(TestParameters(4, 4, 1.0))
        self.parameters.append(TestParameters(19, 19, 1.0))
        self.parameters.append(TestParameters(4, 2, 1.0))
        self.parameters.append(TestParameters(2, 4, 1.0))
        self.parameters.append(TestParameters(4, 4, 0.2))
        self.parameters.append(TestParameters(8, 8, 1.0))

    ##########################################################################
    # Test our ability to read and write matrices.
    #  @param[in] self pointer.
    def test_read(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity, format="csr")
            scipy.io.mmwrite(self.scratch_dir + "/matrix1.mtx",
                             scipy.sparse.csr_matrix(matrix1))
            matrix2 = nt.SparseMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = scipy.io.mmread(self.scratch_dir + "/matrix2.mtx")
            norm = abs(scipy.sparse.linalg.norm(matrix1 - ResultMat))
            self.assertLessEqual(norm, THRESHOLD)
    ##########################################################################
    # Test our ability to read and write matrices.
    #  @param[in] self pointer.
    def test_readsymmetric(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.rows,
                                          param.sparsity, format="csr")
            matrix1 = matrix1 + matrix1.T
            scipy.io.mmwrite(self.scratch_dir + "/matrix1.mtx",
                             scipy.sparse.csr_matrix(matrix1))
            matrix2 = nt.SparseMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = scipy.io.mmread(self.scratch_dir + "/matrix2.mtx")
            norm = abs(scipy.sparse.linalg.norm(matrix1 - ResultMat))
            self.assertLessEqual(norm, THRESHOLD)
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
            scipy.io.mmwrite(self.scratch_dir + "/matrix1.mtx",
                             scipy.sparse.csr_matrix(matrix1))
            scipy.io.mmwrite(self.scratch_dir + "/matrix2.mtx",
                             scipy.sparse.csr_matrix(matrix2))
            alpha = random.uniform(1.0, 2.0)
            CheckMat = alpha * matrix1 + matrix2
            matrix1 = nt.SparseMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2 = nt.SparseMatrix(self.scratch_dir + "/matrix2.mtx")
            matrix2.Increment(matrix1, alpha, 0.0)
            matrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = scipy.io.mmread(self.scratch_dir + "/matrix2.mtx")
            norm = abs(scipy.sparse.linalg.norm(CheckMat - ResultMat))
            self.assertLessEqual(norm, THRESHOLD)
    ##########################################################################
    # Test our ability to dot two matrices.
    #  @param[in] self pointer.

    def test_dot(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            matrix2 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            scipy.io.mmwrite(self.scratch_dir + "/matrix1.mtx",
                             scipy.sparse.csr_matrix(matrix1))
            scipy.io.mmwrite(self.scratch_dir + "/matrix2.mtx",
                             scipy.sparse.csr_matrix(matrix2))
            check = numpy.sum(numpy.multiply(
                matrix1.todense(), matrix2.todense()))
            matrix1 = nt.SparseMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2 = nt.SparseMatrix(self.scratch_dir + "/matrix2.mtx")
            result = matrix2.Dot(matrix1)
            norm = result - check
            self.assertLessEqual(norm, THRESHOLD)
    ##########################################################################
    # Test our ability to transpose a matrix.
    #  @param[in] self pointer.

    def test_transpose(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            scipy.io.mmwrite(self.scratch_dir + "/matrix1.mtx",
                             scipy.sparse.csr_matrix(matrix1))

            matrix2 = nt.SparseMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2T = nt.SparseMatrix(matrix2.GetRows(), matrix2.GetColumns())
            matrix2T.TransposeSparseMatrix(matrix2)
            matrix2T.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")

            CheckMat = matrix1.T
            ResultMat = scipy.io.mmread(self.scratch_dir + "/matrix2.mtx")
            norm = abs(scipy.sparse.linalg.norm(CheckMat - ResultMat))
            self.assertLessEqual(norm, THRESHOLD)

    ##########################################################################
    # Test our ability to pairwise multiply two matrices.
    #  @param[in] self pointer.
    def test_pairwise(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            matrix2 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            scipy.io.mmwrite(self.scratch_dir + "/matrix1.mtx",
                             scipy.sparse.csr_matrix(matrix1))
            scipy.io.mmwrite(self.scratch_dir + "/matrix2.mtx",
                             scipy.sparse.csr_matrix(matrix2))
            CheckMat = numpy.multiply(matrix1.todense(), matrix2.todense())

            ntmatrix1 = nt.SparseMatrix(self.scratch_dir + "/matrix1.mtx")
            ntmatrix2 = nt.SparseMatrix(self.scratch_dir + "/matrix2.mtx")
            ntmatrix3 = nt.SparseMatrix(ntmatrix1.GetColumns(),
                                        ntmatrix1.GetRows())
            ntmatrix3.PairwiseMultiply(ntmatrix1, ntmatrix2)
            ntmatrix3.WriteToMatrixMarket(self.scratch_dir + "/matrix3.mtx")

            ResultMat = scipy.io.mmread(self.scratch_dir + "/matrix3.mtx")
            norm = abs(scipy.linalg.norm(CheckMat - ResultMat))
            self.assertLessEqual(norm, THRESHOLD)

    ##########################################################################
    # Test our ability to multiply two matrices.
    #  @param[in] self pointer.
    def test_multiply(self):
        for param in self.parameters:
            matrix1 = scipy.sparse.random(param.rows, param.columns,
                                          param.sparsity,
                                          format="csr")
            matrix2 = scipy.sparse.random(param.columns, param.columns,
                                          param.sparsity,
                                          format="csr")
            scipy.io.mmwrite(self.scratch_dir + "/matrix1.mtx",
                             scipy.sparse.csr_matrix(matrix1))
            scipy.io.mmwrite(self.scratch_dir + "/matrix2.mtx",
                             scipy.sparse.csr_matrix(matrix2))
            alpha = random.uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > THRESHOLD:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = nt.SparseMatrix(self.scratch_dir + "/matrix1.mtx")
            ntmatrix2 = nt.SparseMatrix(self.scratch_dir + "/matrix2.mtx")
            ntmatrix3 = nt.SparseMatrix(ntmatrix2.GetColumns(),
                                        ntmatrix1.GetRows())
            memory_pool = nt.MatrixMemoryPool(ntmatrix2.GetColumns(),
                                              ntmatrix1.GetRows())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, False, False, alpha, beta, 0.0,
                           memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.scratch_dir + "/matrix3.mtx")

            ResultMat = scipy.io.mmread(self.scratch_dir + "/matrix3.mtx")
            norm = abs(scipy.sparse.linalg.norm(CheckMat - ResultMat))
            self.assertLessEqual(norm, THRESHOLD)


###############################################################################
if __name__ == '__main__':
    unittest.main()
