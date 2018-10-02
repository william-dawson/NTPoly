'''
@package testSparseMatrix
A test suite for local matrices.
'''
import unittest
import NTPolySwig as nt

import scipy
from numpy import sum, multiply, conj
import scipy.sparse
from random import uniform, randint
from scipy.linalg import eigh
from scipy.sparse import random, csr_matrix
from scipy.sparse.linalg import norm
from numpy import diag
from numpy.linalg import norm as normd

from scipy.io import mmwrite, mmread
import os

from helpers import THRESHOLD
from helpers import result_file


class TestParameters:
    '''An internal class for holding test parameters.'''

    def __init__(self, rows, columns, sparsity):
        '''Default constructor
        @param[in] rows matrix rows.
        @param[in] columns matrix columns.
        @param[in] sparsity matrix sparsity.
        '''
        # Matrix rows.
        self.rows = rows
        # Matrix columns.
        self.columns = columns
        # Matrix sparsity.
        self.sparsity = sparsity

    def create_matrix(self, square=False, complex=False):
        '''
        Function to create a matrix for a given set of parameters.
        '''
        r = self.rows
        c = self.columns
        s = self.sparsity
        if square:
            r = c
        if complex:
            mat = random(r, c, s, format="csr")
            mat += 1j * random(r, c, s, format="csr")
        else:
            mat = random(r, c, s, format="csr")

        return csr_matrix(mat)


class TestLocalMatrix(unittest.TestCase):
    '''A test class for local matrices.'''
    # Parameters for the matrices
    parameters = []
    # Location of the scratch directory.
    scratch_dir = os.environ['SCRATCHDIR']
    SMatrix = nt.Matrix_lsr
    MatrixMemoryPool = nt.MatrixMemoryPool_r
    complex = False

    def setUp(self):
        '''Set up a test.'''
        self.parameters = []
        self.parameters.append(TestParameters(2, 4, 0.0))
        self.parameters.append(TestParameters(8, 8, 0.0))
        self.parameters.append(TestParameters(2, 2, 1.0))
        self.parameters.append(TestParameters(4, 4, 1.0))
        self.parameters.append(TestParameters(19, 19, 1.0))
        self.parameters.append(TestParameters(4, 2, 1.0))
        self.parameters.append(TestParameters(2, 4, 1.0))
        self.parameters.append(TestParameters(4, 4, 0.2))
        self.parameters.append(TestParameters(8, 8, 1.0))

    def test_read(self):
        '''Test routines to read and write matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            matrix2 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(matrix1 - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_readsymmetric(self):
        '''Test routines to read and write matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex, square=True)
            matrix1 = matrix1 + matrix1.H
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            matrix2 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(matrix1 - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_addition(self):
        '''Test routines to add together matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            mmwrite(self.scratch_dir + "/matrix2.mtx", matrix2)
            alpha = uniform(1.0, 2.0)
            CheckMat = alpha * matrix1 + matrix2
            matrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2 = self.SMatrix(self.scratch_dir + "/matrix2.mtx")
            matrix2.Increment(matrix1, alpha, 0.0)
            matrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_addzero(self):
        '''Test routines to add together a matrix and zero.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)

            CheckMat = matrix1
            matrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2 = self.SMatrix(matrix1.GetColumns(),
                                   matrix1.GetRows())
            matrix2.Increment(matrix1, 1.0, 0.0)
            matrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_addzeroreverse(self):
        '''Test routines to add together a matrix and zero.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)

            CheckMat = matrix1
            matrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2 = self.SMatrix(matrix1.GetColumns(), matrix1.GetRows())
            matrix1.Increment(matrix2, 1.0, 0.0)
            matrix1.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")
            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_dot(self):
        '''Test routines to dot two matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            mmwrite(self.scratch_dir + "/matrix2.mtx", matrix2)
            check = sum(multiply(matrix1.todense(), matrix2.todense()))
            matrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2 = self.SMatrix(self.scratch_dir + "/matrix2.mtx")
            result = matrix2.Dot(matrix1)
            normval = result - check
            self.assertLessEqual(normval, THRESHOLD)

    def test_transpose(self):
        '''Test routines to transpose a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)

            matrix2 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2T = self.SMatrix(matrix2.GetRows(), matrix2.GetColumns())
            matrix2T.Transpose(matrix2)
            matrix2T.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")

            CheckMat = matrix1.T
            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_pairwise(self):
        '''Test routines to pairwise multiply two matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            mmwrite(self.scratch_dir + "/matrix2.mtx", matrix2)
            CheckMat = csr_matrix(
                multiply(matrix1.todense(), matrix2.todense()))

            ntmatrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            ntmatrix2 = self.SMatrix(self.scratch_dir + "/matrix2.mtx")
            ntmatrix3 = self.SMatrix(
                ntmatrix1.GetColumns(), ntmatrix1.GetRows())
            ntmatrix3.PairwiseMultiply(ntmatrix1, ntmatrix2)
            ntmatrix3.WriteToMatrixMarket(self.scratch_dir + "/matrix3.mtx")

            ResultMat = mmread(self.scratch_dir + "/matrix3.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_multiply(self):
        '''Test routines to multiply two matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex).H
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            mmwrite(self.scratch_dir + "/matrix2.mtx", matrix2)
            alpha = uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > THRESHOLD:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            ntmatrix2 = self.SMatrix(self.scratch_dir + "/matrix2.mtx")
            ntmatrix3 = self.SMatrix(ntmatrix2.GetColumns(),
                                     ntmatrix1.GetRows())
            memory_pool = self.MatrixMemoryPool(ntmatrix2.GetColumns(),
                                                ntmatrix1.GetRows())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, False, False, alpha, beta,
                           0.0, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.scratch_dir + "/matrix3.mtx")

            ResultMat = mmread(self.scratch_dir + "/matrix3.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_multiply_zero(self):
        '''Test routines to multiply two matrices where one is zero.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = 0 * param.create_matrix(complex=self.complex).H
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            mmwrite(self.scratch_dir + "/matrix2.mtx", matrix2)
            alpha = uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > THRESHOLD:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            ntmatrix2 = self.SMatrix(self.scratch_dir + "/matrix2.mtx")
            ntmatrix3 = self.SMatrix(ntmatrix2.GetColumns(),
                                     ntmatrix1.GetRows())
            memory_pool = self.MatrixMemoryPool(ntmatrix2.GetColumns(),
                                                ntmatrix1.GetRows())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, False, False, alpha, beta,
                           0.0, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.scratch_dir + "/matrix3.mtx")

            ResultMat = mmread(self.scratch_dir + "/matrix3.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_eigendecomposition(self):
        '''Test the dense eigen decomposition'''
        for param in self.parameters:
            matrix1 = param.create_matrix(square=True, complex=self.complex)
            matrix1 = matrix1 + matrix1.H
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            w, vdense = eigh(matrix1.todense())
            CheckV = csr_matrix(vdense)

            ntmatrix = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            V = self.SMatrix(ntmatrix.GetColumns(), ntmatrix.GetColumns())

            ntmatrix.EigenDecomposition(V, THRESHOLD)
            V.WriteToMatrixMarket(self.scratch_dir + "/vmat.mtx")

            ResultV = mmread(self.scratch_dir + "/vmat.mtx")
            CheckD = diag((CheckV.H.dot(matrix1).dot(CheckV)).todense())
            ResultD = diag((ResultV.H.dot(matrix1).dot(ResultV)).todense())
            normvalv = abs(normd(CheckD - ResultD))

            self.assertLessEqual(normvalv, THRESHOLD)

    def test_get_row(self):
        '''Test function that extracts a row from the matrix'''
        for param in self.parameters:
            if param.rows == 0:
                continue
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            row_num = randint(0, param.rows - 1)
            CheckMat = matrix1[row_num, :]

            ntmatrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            ntmatrix2 = self.SMatrix(ntmatrix1.GetColumns(), 1)
            ntmatrix1.ExtractRow(row_num, ntmatrix2)
            ntmatrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")

            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_get_column(self):
        '''Test function that extracts a column from the matrix'''
        for param in self.parameters:
            if param.columns == 0:
                continue
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            column_num = randint(0, param.columns - 1)
            CheckMat = matrix1[:, column_num]

            ntmatrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            ntmatrix2 = self.SMatrix(1, ntmatrix1.GetRows())
            ntmatrix1.ExtractColumn(column_num, ntmatrix2)
            ntmatrix2.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")

            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)


class TestLocalMatrix_c(TestLocalMatrix):
    '''Specialization for complex matrices'''
    SMatrix = nt.Matrix_lsc
    MatrixMemoryPool = nt.MatrixMemoryPool_c
    complex = True

    def test_conjugatetranspose(self):
        '''Test routines to compute the conjugate transpose of a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)

            matrix2 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2T = self.SMatrix(matrix2.GetRows(), matrix2.GetColumns())
            matrix2T.Transpose(matrix2)
            matrix2T.Conjugate()
            matrix2T.WriteToMatrixMarket(self.scratch_dir + "/matrix2.mtx")

            CheckMat = matrix1.H
            ResultMat = mmread(self.scratch_dir + "/matrix2.mtx")
            normval = abs(norm(CheckMat - ResultMat))
            self.assertLessEqual(normval, THRESHOLD)

    def test_dot(self):
        '''Test routines to dot two matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.scratch_dir + "/matrix1.mtx", matrix1)
            mmwrite(self.scratch_dir + "/matrix2.mtx", matrix2)
            check = sum(multiply(conj(matrix1.todense()), matrix2.todense()))
            matrix1 = self.SMatrix(self.scratch_dir + "/matrix1.mtx")
            matrix2 = self.SMatrix(self.scratch_dir + "/matrix2.mtx")
            result = matrix2.Dot(matrix1)
            normval = result - check
            self.assertLessEqual(normval, THRESHOLD)


###############################################################################
if __name__ == '__main__':
    unittest.main()
