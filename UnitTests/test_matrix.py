"""
A test suite for local matrices.
"""
import unittest
import NTPolySwig as nt

from scipy.io import mmwrite, mmread


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
        from scipy.sparse import random, csr_matrix

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
    from os import environ
    from os.path import join
    # Parameters for the matrices
    parameters = []
    # Location of the scratch directory.
    scratch_dir = environ['SCRATCHDIR']
    file1 = join(scratch_dir, "matrix1.mtx")
    file2 = join(scratch_dir, "matrix2.mtx")
    file3 = join(scratch_dir, "matrix3.mtx")
    SMatrix = nt.Matrix_lsr
    MatrixMemoryPool = nt.MatrixMemoryPool_r
    complex = False
    TripletList = nt.TripletList_r
    Triplet = nt.Triplet_r

    def _compare_mat(self, val1, val2):
        from helpers import THRESHOLD
        from scipy.sparse.linalg import norm

        normval = abs(norm(val1 - val2))
        self.assertLessEqual(normval, THRESHOLD)

    def _compare(self, val1, val2):
        from helpers import THRESHOLD
        from scipy.linalg import norm

        normval = abs(norm(val1 - val2))
        self.assertLessEqual(normval, THRESHOLD)

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
            mmwrite(self.file1, matrix1)
            matrix2 = self.SMatrix(self.file1)
            matrix2.WriteToMatrixMarket(self.file2)
            ResultMat = mmread(self.file2)
            self._compare_mat(matrix1, ResultMat)

    def test_readcircular(self):
        '''Test routines to read a matrix produced by ntpoly.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)
            matrix2 = self.SMatrix(self.file1)
            matrix2.WriteToMatrixMarket(self.file2)
            matrix3 = self.SMatrix(self.file2)
            matrix3.WriteToMatrixMarket(self.file3)
            ResultMat = mmread(self.file3)

            self._compare_mat(matrix1, ResultMat)

    def test_readsymmetric(self):
        '''Test routines to read and write matrices.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex, square=True)
            matrix1 = matrix1 + matrix1.H
            mmwrite(self.file1, matrix1)
            matrix2 = self.SMatrix(self.file1)
            matrix2.WriteToMatrixMarket(self.file2)
            ResultMat = mmread(self.file2)

            self._compare_mat(matrix1, ResultMat)

    def test_addition(self):
        '''Test routines to add together matrices.'''
        from random import uniform
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)
            mmwrite(self.file2, matrix2)
            alpha = uniform(1.0, 2.0)
            CheckMat = alpha * matrix1 + matrix2
            matrix1 = self.SMatrix(self.file1)
            matrix2 = self.SMatrix(self.file2)
            matrix2.Increment(matrix1, alpha, 0.0)
            matrix2.WriteToMatrixMarket(self.file2)
            ResultMat = mmread(self.file2)

            self._compare_mat(CheckMat, ResultMat)

    def test_addzero(self):
        '''Test routines to add together a matrix and zero.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)

            CheckMat = matrix1
            matrix1 = self.SMatrix(self.file1)
            matrix2 = self.SMatrix(matrix1.GetColumns(), matrix1.GetRows())
            matrix2.Increment(matrix1, 1.0, 0.0)
            matrix2.WriteToMatrixMarket(self.file2)
            ResultMat = mmread(self.file2)
            self._compare_mat(CheckMat, ResultMat)

    def test_addzeroreverse(self):
        '''Test routines to add together a matrix and zero.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)

            CheckMat = matrix1
            matrix1 = self.SMatrix(self.file1)
            matrix2 = self.SMatrix(matrix1.GetColumns(), matrix1.GetRows())
            matrix1.Increment(matrix2, 1.0, 0.0)
            matrix1.WriteToMatrixMarket(self.file2)
            ResultMat = mmread(self.file2)
            self._compare_mat(CheckMat, ResultMat)

    def test_dot(self):
        '''Test routines to dot two matrices.'''
        from numpy import sum, multiply
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)
            mmwrite(self.file2, matrix2)
            check = sum(multiply(matrix1.todense(), matrix2.todense()))
            matrix1 = self.SMatrix(self.file1)
            matrix2 = self.SMatrix(self.file2)
            result = matrix2.Dot(matrix1)

            self._compare(result, check)

    def test_transpose(self):
        '''Test routines to transpose a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)

            matrix2 = self.SMatrix(self.file1)
            matrix2T = self.SMatrix(matrix2.GetRows(), matrix2.GetColumns())
            matrix2T.Transpose(matrix2)
            matrix2T.WriteToMatrixMarket(self.file2)

            CheckMat = matrix1.T
            ResultMat = mmread(self.file2)

            self._compare_mat(CheckMat, ResultMat)

    def test_pairwise(self):
        '''Test routines to pairwise multiply two matrices.'''
        from scipy.sparse import csr_matrix
        from numpy import multiply
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)
            mmwrite(self.file2, matrix2)
            CheckMat = csr_matrix(
                multiply(matrix1.todense(), matrix2.todense()))

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(self.file2)
            ntmatrix3 = self.SMatrix(
                ntmatrix1.GetColumns(), ntmatrix1.GetRows())
            ntmatrix3.PairwiseMultiply(ntmatrix1, ntmatrix2)
            ntmatrix3.WriteToMatrixMarket(self.file3)

            ResultMat = mmread(self.file3)
            self._compare_mat(CheckMat, ResultMat)

    def test_multiply(self):
        '''Test routines to multiply two matrices.'''
        from random import uniform
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex).H
            mmwrite(self.file1, matrix1)
            mmwrite(self.file2, matrix2)
            alpha = uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > 0.0:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(self.file2)
            ntmatrix3 = self.SMatrix(ntmatrix2.GetColumns(),
                                     ntmatrix1.GetRows())
            memory_pool = self.MatrixMemoryPool(ntmatrix2.GetColumns(),
                                                ntmatrix1.GetRows())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, False, False, alpha, beta,
                           0.0, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.file3)

            ResultMat = mmread(self.file3)
            self._compare_mat(CheckMat, ResultMat)

    def test_multiply_nt(self):
        '''Test routines to multiply two matrices.'''
        from random import uniform
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex).H
            mmwrite(self.file1, matrix1)
            mmwrite(self.file2, matrix2.T)
            alpha = uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > 0.0:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(self.file2)
            ntmatrix3 = self.SMatrix(ntmatrix2.GetRows(),
                                     ntmatrix1.GetRows())
            memory_pool = self.MatrixMemoryPool(ntmatrix2.GetRows(),
                                                ntmatrix1.GetRows())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, False, True, alpha, beta,
                           0.0, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.file3)

            ResultMat = mmread(self.file3)
            self._compare_mat(CheckMat, ResultMat)

    def test_multiply_tn(self):
        '''Test routines to multiply two matrices.'''
        from random import uniform
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex).H
            mmwrite(self.file1, matrix1.T)
            mmwrite(self.file2, matrix2)
            alpha = uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > 0.0:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(self.file2)
            ntmatrix3 = self.SMatrix(ntmatrix2.GetColumns(),
                                     ntmatrix1.GetColumns())
            memory_pool = self.MatrixMemoryPool(ntmatrix2.GetColumns(),
                                                ntmatrix1.GetColumns())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, True, False, alpha, beta,
                           0.0, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.file3)

            ResultMat = mmread(self.file3)
            self._compare_mat(CheckMat, ResultMat)

    def test_multiply_tt(self):
        '''Test routines to multiply two matrices.'''
        from random import uniform
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex).H
            mmwrite(self.file1, matrix1.T)
            mmwrite(self.file2, matrix2.T)
            alpha = uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > 0.0:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(self.file2)
            ntmatrix3 = self.SMatrix(ntmatrix2.GetRows(),
                                     ntmatrix1.GetColumns())
            memory_pool = self.MatrixMemoryPool(ntmatrix2.GetRows(),
                                                ntmatrix1.GetColumns())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, True, True, alpha, beta,
                           0.0, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.file3)

            ResultMat = mmread(self.file3)
            self._compare_mat(CheckMat, ResultMat)

    def test_multiply_zero(self):
        '''Test routines to multiply two matrices where one is zero.'''
        from random import uniform
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = 0 * param.create_matrix(complex=self.complex).H
            mmwrite(self.file1, matrix1)
            mmwrite(self.file2, matrix2)
            alpha = uniform(1.0, 2.0)
            beta = 0.0
            if abs(beta) > 0.0:
                CheckMat = alpha * matrix1.dot(matrix2) + beta * matrix1
            else:
                CheckMat = alpha * matrix1.dot(matrix2)

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(self.file2)
            ntmatrix3 = self.SMatrix(ntmatrix2.GetColumns(),
                                     ntmatrix1.GetRows())
            memory_pool = self.MatrixMemoryPool(ntmatrix2.GetColumns(),
                                                ntmatrix1.GetRows())
            ntmatrix3.Gemm(ntmatrix1, ntmatrix2, False, False, alpha, beta,
                           0.0, memory_pool)
            ntmatrix3.WriteToMatrixMarket(self.file3)

            ResultMat = mmread(self.file3)
            self._compare_mat(CheckMat, ResultMat)

    def test_get_row(self):
        '''Test function that extracts a row from the matrix'''
        from random import randint
        for param in self.parameters:
            if param.rows == 0:
                continue
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)
            row_num = randint(0, param.rows - 1)
            CheckMat = matrix1[row_num, :]

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(ntmatrix1.GetColumns(), 1)
            ntmatrix1.ExtractRow(row_num, ntmatrix2)
            ntmatrix2.WriteToMatrixMarket(self.file2)

            ResultMat = mmread(self.file2)
            self._compare_mat(CheckMat, ResultMat)

    def test_get_column(self):
        '''Test function that extracts a column from the matrix'''
        from random import randint
        for param in self.parameters:
            if param.columns == 0:
                continue
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)
            column_num = randint(0, param.columns - 1)
            CheckMat = matrix1[:, column_num]

            ntmatrix1 = self.SMatrix(self.file1)
            ntmatrix2 = self.SMatrix(1, ntmatrix1.GetRows())
            ntmatrix1.ExtractColumn(column_num, ntmatrix2)
            ntmatrix2.WriteToMatrixMarket(self.file2)

            ResultMat = mmread(self.file2)
            self._compare_mat(CheckMat, ResultMat)

    def test_scalediag(self):
        '''Test routines to scale by a diagonal matrix.'''
        from copy import deepcopy
        for param in self.parameters:
            matrix = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix)
            CheckMat = deepcopy(matrix)

            tlist = self.TripletList()
            for i in range(matrix.shape[1]):
                t = self.Triplet()
                t.index_column = i + 1
                t.index_row = i + 1
                t.point_value = i
                tlist.Append(t)
                CheckMat[:, i] *= i
            ntmatrix = self.SMatrix(self.file1)
            ntmatrix.DiagonalScale(tlist, 0)
            ntmatrix.WriteToMatrixMarket(self.file2)

            ResultMat = mmread(self.file2)
            self._compare_mat(CheckMat, ResultMat)


class TestLocalMatrix_c(TestLocalMatrix):
    '''Specialization for complex matrices'''
    SMatrix = nt.Matrix_lsc
    MatrixMemoryPool = nt.MatrixMemoryPool_c
    complex = True
    TripletList = nt.TripletList_c
    Triplet = nt.Triplet_c

    def test_conjugatetranspose(self):
        '''Test routines to compute the conjugate transpose of a matrix.'''
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)

            matrix2 = self.SMatrix(self.file1)
            matrix2T = self.SMatrix(matrix2.GetRows(), matrix2.GetColumns())
            matrix2T.Transpose(matrix2)
            matrix2T.Conjugate()
            matrix2T.WriteToMatrixMarket(self.file2)

            CheckMat = matrix1.H
            ResultMat = mmread(self.file2)
            self._compare_mat(CheckMat, ResultMat)

    def test_dot(self):
        '''Test routines to dot two matrices.'''
        from numpy import sum, multiply, conj
        for param in self.parameters:
            matrix1 = param.create_matrix(complex=self.complex)
            matrix2 = param.create_matrix(complex=self.complex)
            mmwrite(self.file1, matrix1)
            mmwrite(self.file2, matrix2)
            check = sum(multiply(conj(matrix1.todense()), matrix2.todense()))
            matrix1 = self.SMatrix(self.file1)
            matrix2 = self.SMatrix(self.file2)
            result = matrix1.Dot(matrix2)

            self._compare(result, check)


###############################################################################
if __name__ == '__main__':
    unittest.main()
