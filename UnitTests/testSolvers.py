##########################################################################
# @package testDistributedSparseMatrix
#  A test suite for the Distributed Sparse Matrix module.
import unittest
import NTPolySwig as nt

import scipy
import scipy.sparse
import scipy.io
import numpy
import random
import time
import os
import math
import numpy.polynomial.chebyshev
from mpi4py import MPI
comm = MPI.COMM_WORLD

from Helpers import THRESHOLD
from Helpers import result_file
from Helpers import scratch_dir

##########################################################################
# A test class for the distributed matrix module.


class TestSolvers(unittest.TestCase):
    # Parameters for the tests
    input_file = scratch_dir + "/input.mtx"
    input_file2 = scratch_dir + "/input2.mtx"
    CheckMat = 0
    my_rank = 0
    matrix_dimension = 32

    ##########################################################################
    # set up all of the tests
    #  @param[in] self pointer.
    @classmethod
    def setUpClass(self):
        rows = int(os.environ['PROCESS_ROWS'])
        columns = int(os.environ['PROCESS_COLUMNS'])
        slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(rows, columns, slices)

    ##########################################################################
    # set up an individual test
    #  @param[in] self pointer
    def setUp(self):
        self.my_rank = comm.Get_rank()
        self.iterative_solver_parameters = nt.IterativeSolverParameters()
        self.fixed_solver_parameters = nt.FixedSolverParameters()
        self.fixed_solver_parameters.SetVerbosity(True)
        self.iterative_solver_parameters.SetVerbosity(True)

    ##########################################################################
    def check_result(self):
        norm = 0
        if (self.my_rank == 0):
            ResultMat = scipy.io.mmread(result_file)
            norm = abs(scipy.sparse.linalg.norm(self.CheckMat - ResultMat))
            print "Norm:", norm
        global_norm = comm.bcast(norm, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    ##########################################################################
    # Test our ability to invert matrices.
    #  @param[in] self pointer.
    def test_invert(self):
        # Starting Matrix
        temp_mat = numpy.zeros((self.matrix_dimension, self.matrix_dimension))
        for j in range(0, self.matrix_dimension):
            for i in range(0, self.matrix_dimension):
                if i == j:
                    temp_mat[j, i] = 1.0
                else:
                    temp_mat[j, i] = random.random() / (float(i - j)**2)
        temp_mat = (temp_mat.T + temp_mat) * 0.5
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        self.CheckMat = scipy.sparse.linalg.inv(
            scipy.sparse.csc_matrix(matrix1))
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        overlap_matrix = nt.DistributedSparseMatrix(
            self.input_file, False)
        inverse_matrix = nt.DistributedSparseMatrix(
            overlap_matrix.GetActualDimension())
        permutation = nt.Permutation(overlap_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)

        nt.InverseSolvers.Invert(overlap_matrix, inverse_matrix,
                                 self.iterative_solver_parameters)

        inverse_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute the inverse square root of matrices.
    #  @param[in] self pointer.
    def test_inversesquareroot(self):
        # Starting Matrix. Care taken to make sure eigenvalues are positive.
        temp_mat = scipy.sparse.rand(self.matrix_dimension, self.matrix_dimension,
                                     density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        identity_matrix = scipy.sparse.identity(self.matrix_dimension)
        temp_mat = temp_mat + identity_matrix * self.matrix_dimension
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: 1.0/numpy.sqrt(x))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        overlap_matrix = nt.DistributedSparseMatrix(
            self.input_file, False)
        inverse_matrix = nt.DistributedSparseMatrix(
            overlap_matrix.GetActualDimension())
        permutation = nt.Permutation(overlap_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix, inverse_matrix,
                                               self.iterative_solver_parameters)
        inverse_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute the square root of matrices.
    #  @param[in] self pointer.
    def test_squareroot(self):
        # Starting Matrix. Care taken to make sure eigenvalues are positive.
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        identity_matrix = scipy.sparse.identity(self.matrix_dimension)
        temp_mat = temp_mat + identity_matrix * self.matrix_dimension
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.sqrt(x))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        root_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.SquareRootSolvers.SquareRoot(input_matrix, root_matrix,
                                        self.iterative_solver_parameters)
        root_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute  general matrix inverse root.
    #  @param[in] self pointer.
    def test_inverseroot(self):
        roots = [1, 2, 4, 5, 6, 7, 8]
        for root in roots:
            print "Root:", root
            # Starting Matrix. Care taken to make sure eigenvalues are
            # positive.
            temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                         self.matrix_dimension,
                                         density=1.0)
            temp_mat = (temp_mat.T + temp_mat)
            identity_matrix = scipy.sparse.identity(self.matrix_dimension)
            temp_mat = temp_mat + identity_matrix * self.matrix_dimension
            matrix1 = scipy.sparse.csr_matrix(temp_mat)

            # Check Matrix
            dense_check = scipy.linalg.funm(temp_mat.todense(),
                                            lambda x: numpy.power(x,-1.0/root))
            self.CheckMat = scipy.sparse.csr_matrix(dense_check)
            if self.my_rank == 0:
                scipy.io.mmwrite(self.input_file,
                                 scipy.sparse.csr_matrix(matrix1))
            comm.barrier()

            # Result Matrix
            input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
            inverse_matrix = nt.DistributedSparseMatrix(
                input_matrix.GetActualDimension())
            permutation = nt.Permutation(input_matrix.GetLogicalDimension())
            permutation.SetRandomPermutation()
            self.iterative_solver_parameters.SetLoadBalance(permutation)
            nt.RootSolvers.ComputeInverseRoot(input_matrix, inverse_matrix,
                                              root,
                                              self.iterative_solver_parameters)
            inverse_matrix.WriteToMatrixMarket(result_file)
            comm.barrier()

            self.check_result()

    ##########################################################################
    # Test our ability to compute  general matrix root.
    #  @param[in] self pointer.
    def test_root(self):
        roots = [1, 2, 4, 5, 6, 7, 8]
        for root in roots:
            print "Root", root
            # Starting Matrix. Care taken to make sure eigenvalues are
            # positive.
            temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                         self.matrix_dimension,
                                         density=1.0)
            temp_mat = (temp_mat.T + temp_mat)
            identity_matrix = scipy.sparse.identity(self.matrix_dimension)
            temp_mat = temp_mat + identity_matrix * self.matrix_dimension
            matrix1 = scipy.sparse.csr_matrix(temp_mat)

            # Check Matrix
            dense_check = scipy.linalg.funm(temp_mat.todense(),
                                            lambda x: numpy.power(x,1.0/root))
            self.CheckMat = scipy.sparse.csr_matrix(dense_check)
            if self.my_rank == 0:
                scipy.io.mmwrite(self.input_file,
                                 scipy.sparse.csr_matrix(matrix1))
            comm.barrier()

            # Result Matrix
            input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
            inverse_matrix = nt.DistributedSparseMatrix(
                input_matrix.GetActualDimension())
            permutation = nt.Permutation(input_matrix.GetLogicalDimension())
            permutation.SetRandomPermutation()
            self.iterative_solver_parameters.SetLoadBalance(permutation)
            nt.RootSolvers.ComputeRoot(input_matrix, inverse_matrix,
                                       root,
                                       self.iterative_solver_parameters)
            inverse_matrix.WriteToMatrixMarket(result_file)
            comm.barrier()
            self.check_result()

    ##########################################################################
    # Test our ability to compute the matrix sign function.
    #  @param[in] self pointer.
    def test_signfunction(self):
        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.sign(x))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        sign_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.SignSolvers.Compute(input_matrix, sign_matrix,
                               self.iterative_solver_parameters)
        sign_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute the matrix exponential.
    #  @param[in] self pointer.
    def test_exponentialfunction(self):
        # Starting Matrix
        temp_mat = scipy.sparse.rand(
            self.matrix_dimension, self.matrix_dimension, density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.exp(x))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        exp_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.ExponentialSolvers.ComputeExponential(input_matrix, exp_matrix,
                                                 self.fixed_solver_parameters)
        exp_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute the matrix logarithm.
    #  @param[in] self pointer.
    def test_logarithmfunction(self):
        # Starting Matrix. Care taken to make sure eigenvalues are positive.
        temp_mat = scipy.sparse.rand(self.matrix_dimension, self.matrix_dimension,
                                     density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        identity_matrix = scipy.sparse.identity(self.matrix_dimension)
        temp_mat = temp_mat + identity_matrix * self.matrix_dimension
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.log(x))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        log_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.ExponentialSolvers.ComputeLogarithm(input_matrix, log_matrix,
                                               self.fixed_solver_parameters)
        log_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute the matrix sine.
    #  @param[in] self pointer.
    def test_sinfunction(self):
        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * \
            (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.sin(x))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        sin_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.TrigonometrySolvers.Sine(input_matrix, sin_matrix,
                                    self.fixed_solver_parameters)
        sin_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute the matrix cosine.
    #  @param[in] self pointer.
    def test_cosfunction(self):
        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension, density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * \
            (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.cos(x))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        cos_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.TrigonometrySolvers.Cosine(input_matrix, cos_matrix,
                                      self.fixed_solver_parameters)
        cos_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute a matrix polynomial using horner's method.
    #  @param[in] self pointer.
    def test_hornerfunction(self):
        # Coefficients of the polynomial
        coef = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]

        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        val, vec = numpy.linalg.eigh(temp_mat.todense())
        for i in range(0, len(val)):
            temp = val[i]
            val[i] = 0
            for j in range(0, len(coef)):
                val[i] = val[i] + coef[j] * (temp**j)
        temp_poly = numpy.dot(numpy.dot(vec, numpy.diag(val)), vec.T)
        self.CheckMat = scipy.sparse.csr_matrix(temp_poly)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.Polynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.HornerCompute(input_matrix, poly_matrix,
                                 self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute a matrix polynomial using the paterson
    #  stockmeyer method.
    #  @param[in] self pointer.
    def test_patersonstockmeyerfunction(self):
        # Coefficients of the polynomial
        coef = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]

        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Check Matrix
        val, vec = numpy.linalg.eigh(temp_mat.todense())
        for i in range(0, len(val)):
            temp = val[i]
            val[i] = 0
            for j in range(0, len(coef)):
                val[i] = val[i] + coef[j] * (temp**j)
        temp_poly = numpy.dot(numpy.dot(vec, numpy.diag(val)), vec.T)
        self.CheckMat = scipy.sparse.csr_matrix(temp_poly)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.Polynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.PatersonStockmeyerCompute(input_matrix, poly_matrix,
                                             self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute using Chebyshev polynomials.
    #  @param[in] self pointer.
    def test_chebyshevfunction(self):
        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Function
        x = numpy.linspace(-1.0, 1.0, 200)
        y = [math.cos(i) + math.sin(i) for i in x]
        coef = numpy.polynomial.chebyshev.chebfit(x, y, 10)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.polynomial.chebyshev.chebval(x, coef))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.ChebyshevPolynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.Compute(input_matrix, poly_matrix,
                           self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to compute using Chebyshev polynomials recursively.
    #  @param[in] self pointer.
    def test_recursivechebyshevfunction(self):
        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = scipy.sparse.csr_matrix(temp_mat)

        # Function
        x = numpy.linspace(-1.0, 1.0, 200)
        # y = [scipy.special.erfc(i) for i in x]
        y = [math.exp(i) for i in x]
        coef = numpy.polynomial.chebyshev.chebfit(x, y, 16 - 1)

        # Check Matrix
        dense_check = scipy.linalg.funm(temp_mat.todense(),
                                        lambda x: numpy.polynomial.chebyshev.chebval(x, coef))
        self.CheckMat = scipy.sparse.csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file,
                             scipy.sparse.csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.ChebyshevPolynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.ComputeFactorized(
            input_matrix, poly_matrix, self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    ##########################################################################
    # Test our ability to solve general matrix equations.
    #  @param[in] self pointer.
    def test_cgsolve(self):
        # Starting Matrix
        A = scipy.sparse.rand(self.matrix_dimension, self.matrix_dimension,
                              density=1.0)
        B = scipy.sparse.rand(self.matrix_dimension, self.matrix_dimension,
                              density=1.0)
        A = A.T.dot(A)

        # Check Matrix
        Ainv = scipy.sparse.linalg.inv(A)
        self.CheckMat = scipy.sparse.csr_matrix(Ainv.dot(B))
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file, scipy.sparse.csr_matrix(A))
            scipy.io.mmwrite(self.input_file2, scipy.sparse.csr_matrix(B))
        comm.barrier()

        # Result Matrix
        AMat = nt.DistributedSparseMatrix(self.input_file, False)
        XMat = nt.DistributedSparseMatrix(AMat.GetActualDimension())
        BMat = nt.DistributedSparseMatrix(self.input_file2, False)
        permutation = nt.Permutation(AMat.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)

        nt.LinearSolvers.CGSolver(AMat, XMat, BMat,
                                  self.iterative_solver_parameters)

        XMat.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()


if __name__ == '__main__':
    unittest.main()
