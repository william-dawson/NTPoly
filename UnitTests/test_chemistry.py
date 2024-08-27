"""
A test suite for chemistry focused solvers.
"""
import unittest
import NTPolySwig as nt
from mpi4py import MPI
# MPI global communicator
comm = MPI.COMM_WORLD


class TestChemistry:
    '''A test class for the chemistry solvers.'''
    # Matrix to compare to
    CheckMat = 0
    # Rank of the current process
    my_rank = 0

    def create_matrices(self):
        '''
        Create the test matrix with the following parameters.
        '''
        from numpy import diag, sqrt
        from scipy.linalg import eigh, funm
        from scipy.sparse import csr_matrix, rand

        fock = rand(self.mat_dim, self.mat_dim, density=1.0)
        if self.is_complex:
            fock += 1j * rand(self.mat_dim, self.mat_dim, density=1.0)
            fock = fock + fock.getH()
        else:
            fock = fock + fock.T
        overlap = rand(self.mat_dim, self.mat_dim, density=1.0)
        overlap = overlap.T.dot(overlap)

        # Make sure the overlap is well conditioned.
        w, v = eigh(overlap.todense())
        w += 0.2
        overlap = csr_matrix(v.T.dot(diag(w).dot(v)))

        isq = funm(overlap.todense(), lambda x: 1.0 / sqrt(x))
        wfock = isq.dot(fock.todense()).dot(isq)

        # Add a gap
        w, v = eigh(wfock)
        gap = (w[-1] - w[0]) / 2.0
        w[int(self.nel):] += gap
        if self.is_complex:
            wfock = v.conj().T.dot(diag(w).dot(v))
        else:
            wfock = v.T.dot(diag(w).dot(v))

        # Compute the density
        w[:int(self.nel)] = 2.0
        w[int(self.nel):] = 0.0
        if self.is_complex:
            density = isq.dot(v.dot(diag(w).dot(v.conj().T))).dot(isq)
        else:
            density = isq.dot(v.dot(diag(w).dot(v.T))).dot(isq)

        self.write_matrix(fock, self.hamiltonian)
        self.write_matrix(overlap, self.overlap)
        self.write_matrix(density, self.density)

    def write_matrix(self, mat, file_name):
        from scipy.io import mmwrite
        from scipy.sparse import csr_matrix

        if self.my_rank == 0:
            mmwrite(file_name, csr_matrix(mat))
        comm.barrier()

    @classmethod
    def setUpClass(self):
        '''Set up all of the tests.'''
        from os import environ
        rows = int(environ['PROCESS_ROWS'])
        columns = int(environ['PROCESS_COLUMNS'])
        slices = int(environ['PROCESS_SLICES'])
        nt.ConstructGlobalProcessGrid(rows, columns, slices)

    @classmethod
    def tearDownClass(self):
        '''Cleanup this test'''
        nt.DestructGlobalProcessGrid()

    def tearDown(self):
        from helpers import log_file
        from yaml import load, dump, SafeLoader
        from sys import stdout
        if nt.GetGlobalIsRoot():
            nt.DeactivateLogger()
            with open(log_file) as ifile:
                data = load(ifile, Loader=SafeLoader)
            dump(data, stdout)

    def setUp(self):
        '''Set up an individual test.'''
        from os import environ
        from os.path import join
        from helpers import scratch_dir, log_file

        self.my_rank = comm.Get_rank()
        self.solver_parameters = nt.SolverParameters()
        self.solver_parameters.SetVerbosity(True)
        self.geomh1 = environ["GEOMH1"]
        self.geomo1 = environ["GEOMO1"]
        self.geomo2 = environ["GEOMO2"]
        self.geomd2 = environ["GEOMD2"]
        self.realio = environ["REALIO"]
        self.nel = 5.0

        self.hamiltonian = join(scratch_dir, "rf.mtx")
        self.overlap = join(scratch_dir, "rs.mtx")
        self.density = join(scratch_dir, "rd.mtx")
        self.mat_dim = 7

        if nt.GetGlobalIsRoot():
            nt.ActivateLogger(log_file, True)

    def check_full(self):
        '''Compare two computed matrices.'''
        from helpers import THRESHOLD, result_file
        from scipy.sparse.linalg import norm
        from scipy.io import mmread

        normval = 0
        if (self.my_rank == 0):
            ResultMat = 2.0 * mmread(result_file)
            self.CheckMat = mmread(self.density)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)

    def check_full_extrap(self):
        '''Compare two computed matrices.'''
        from helpers import EXTRAPTHRESHOLD, result_file
        from scipy.sparse.linalg import norm
        from scipy.io import mmread

        normval = 0
        if (self.my_rank == 0):
            ResultMat = 2.0 * mmread(result_file)
            self.CheckMat = mmread(self.geomd2)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, EXTRAPTHRESHOLD)

    def compute_cp(self):
        '''Compute the chemical potential, homo, lumo'''
        from scipy.io import mmread
        from scipy.linalg import eigh

        fock_matrix = mmread(self.hamiltonian)
        overlap_matrix = mmread(self.overlap)
        eig_vals = eigh(a=fock_matrix.todense(), b=overlap_matrix.todense(),
                        eigvals_only=True)
        homo = int(self.nel) - 1
        lumo = homo + 1
        cp = eig_vals[homo] + (eig_vals[lumo] - eig_vals[homo]) / 2.0
        return cp, eig_vals[homo], eig_vals[lumo]

    def check_energy(self, energy):
        '''Compute the chemical potential, homo, lumo'''
        from scipy.io import mmread
        from helpers import THRESHOLD
        from scipy.linalg import eigh
        from numpy import floor, ceil

        fock_matrix = mmread(self.hamiltonian)
        overlap_matrix = mmread(self.overlap)

        eig_vals = eigh(a=fock_matrix.todense(), b=overlap_matrix.todense(),
                        eigvals_only=True)

        computed = 0
        for i, v in enumerate(eig_vals):
            if i < floor(self.nel):
                computed += v
            elif ceil(self.nel) == i:
                computed += (ceil(self.nel) - self.nel) * v

        self.assertLessEqual(abs(energy - computed), THRESHOLD)

    def check_cp(self, computed):
        '''Compare two computed chemical potentials.'''
        cp, homo, lumo = self.compute_cp()

        if computed > homo and computed < lumo:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def basic_solver(self, SRoutine, cpcheck=True, temp=None):
        '''Test various kinds of density matrix solvers.'''
        from helpers import result_file

        self.create_matrices()
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        overlap_matrix = nt.Matrix_ps(self.overlap)
        inverse_sqrt_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        density_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        if temp is None:
            energy_value, chemical_potential = SRoutine(fock_matrix,
                                                        inverse_sqrt_matrix,
                                                        self.nel,
                                                        density_matrix,
                                                        self.solver_parameters)
        else:
            inv_temp = 1/(temp * 3.166811563*10**(-6))
            print("::: Temperature", temp, inv_temp, self.nel)
            energy_value, chemical_potential = SRoutine(fock_matrix,
                                                        inverse_sqrt_matrix,
                                                        self.nel,
                                                        density_matrix,
                                                        inv_temp,
                                                        self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full()
        if cpcheck:
            self.check_cp(chemical_potential)
        self.check_energy(energy_value)
        comm.barrier()

    def test_scaleandfold(self):
        '''Test the scale and fold method.'''
        from helpers import result_file

        SRoutine = nt.DensityMatrixSolvers.ScaleAndFold
        self.create_matrices()
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        overlap_matrix = nt.Matrix_ps(self.overlap)
        inverse_sqrt_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        density_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)

        # This method needs the homo and lumo
        cp, homo, lumo = self.compute_cp()
        print(":::::", homo, lumo, cp)
        SRoutine(fock_matrix, inverse_sqrt_matrix, self.nel, density_matrix,
                 homo, lumo, self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full()
        comm.barrier()

    def test_pm(self):
        '''Test our ability to compute the density matrix with PM.'''
        self.basic_solver(nt.DensityMatrixSolvers.PM)

    def test_trs2(self):
        '''Test routines to compute the density matrix with TRS2.'''
        self.basic_solver(nt.DensityMatrixSolvers.TRS2)

    def test_trs4(self):
        '''Test routines to compute the density matrix with TRS4.'''
        self.basic_solver(nt.DensityMatrixSolvers.TRS4)

    def test_HPCP(self):
        '''Test routines to compute the density matrix with HPCP.'''
        self.basic_solver(nt.DensityMatrixSolvers.HPCP)

    def test_densedensity(self):
        '''Test routines to compute the density matrix with Dense Method.'''
        self.basic_solver(nt.DensityMatrixSolvers.DenseDensity)

    def test_foe_low(self):
        '''Test the fermi operator expansion at a low temperature.'''
        self.basic_solver(nt.FermiOperator.ComputeDenseFOE, temp=50)

    def test_wom_gc(self):
        '''Test the wom_gc solver.'''
        from helpers import result_file, THRESHOLD
        from scipy.io import mmread
        from scipy.linalg import eigh, funm
        from numpy.linalg import norm
        from numpy import exp, diag, sqrt, trace, matrix

        SRoutine = nt.FermiOperator.WOM_GC
        beta = 105

        self.create_matrices()
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        overlap_matrix = nt.Matrix_ps(self.overlap)
        inverse_sqrt_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        density_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)
        self.solver_parameters.SetStepThreshold(1e-5)

        # Reference Solution
        F = mmread(self.hamiltonian).todense()
        S = mmread(self.overlap).todense()
        ISQ = funm(S, lambda x: 1/sqrt(x))
        FOrth = ISQ.dot(F).dot(ISQ)

        w, v = eigh(FOrth)
        mu = w[int(self.nel)-1] + 0.5*(w[int(self.nel)] - w[int(self.nel)-1])
        we = [1.0/(1 + exp(beta*(x - mu))) for x in w]

        if self.is_complex:
            DOrth = v @ diag(we) @ matrix(v).getH()
        else:
            DOrth = v @ diag(we) @ v.T
        D = ISQ.dot(DOrth).dot(ISQ)
        check_energy = trace(FOrth.dot(DOrth))

        # Solve with NTPoly
        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        energy_value = SRoutine(fock_matrix, inverse_sqrt_matrix,
                                density_matrix, mu, beta,
                                self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)

        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            self.CheckMat = D
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        global_en = comm.bcast(check_energy - energy_value, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        self.assertLessEqual(global_en, THRESHOLD)

    def test_wom_c(self):
        '''Test the wom_c solver.'''
        from helpers import result_file, THRESHOLD
        from scipy.io import mmread
        from scipy.linalg import eigh, funm
        from numpy.linalg import norm
        from numpy import exp, diag, sqrt, trace, matrix

        SRoutine = nt.FermiOperator.WOM_C
        beta = 105

        self.create_matrices()
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        overlap_matrix = nt.Matrix_ps(self.overlap)
        inverse_sqrt_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        density_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())

        permutation = nt.Permutation(fock_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)
        self.solver_parameters.SetStepThreshold(1e-5)

        # Reference Solution
        F = mmread(self.hamiltonian).todense()
        S = mmread(self.overlap).todense()
        ISQ = funm(S, lambda x: 1/sqrt(x))
        FOrth = ISQ.dot(F).dot(ISQ)

        w, v = eigh(FOrth)
        mu = w[int(self.nel)-1] + 0.5*(w[int(self.nel)] - w[int(self.nel)-1])
        we = [1.0/(1 + exp(beta*(x - mu))) for x in w]
        DOrth = 1 * S
        while abs(trace(DOrth) - self.nel) > 1e-12:
            we = [1.0/(1 + exp(beta*(x - mu))) for x in w]
            if self.is_complex:
                DOrth = v @ diag(we) @ matrix(v).getH()
            else:
                DOrth = v @ diag(we) @ v.T
        D = ISQ.dot(DOrth).dot(ISQ)
        check_energy = trace(FOrth.dot(DOrth))

        # Solve with NTPoly
        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix,
                                               inverse_sqrt_matrix,
                                               self.solver_parameters)
        energy_value = SRoutine(fock_matrix, inverse_sqrt_matrix,
                                density_matrix, self.nel, beta,
                                self.solver_parameters)

        density_matrix.WriteToMatrixMarket(result_file)

        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            self.CheckMat = D
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        global_en = comm.bcast(check_energy - energy_value, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        self.assertLessEqual(global_en, THRESHOLD)

    def test_mcweeny_step(self):
        from scipy.io import mmread
        from helpers import result_file, THRESHOLD
        from scipy.sparse.linalg import norm

        # Reference Solution
        dmat = mmread(self.density)
        smat = mmread(self.overlap)
        result = 3 * dmat.dot(dmat) - 2 * dmat.dot(dmat).dot(dmat)
        result_s = 3 * dmat.dot(smat).dot(dmat) - \
            2 * dmat.dot(smat).dot(dmat).dot(smat).dot(dmat)

        # NTPoly
        d_matrix = nt.Matrix_ps(self.density)
        dout_matrix = nt.Matrix_ps(d_matrix.GetActualDimension())
        nt.DensityMatrixSolvers.McWeenyStep(d_matrix, dout_matrix)
        dout_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        # Compare
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            normval = abs(norm(result - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        comm.barrier()

        # Overlap Version
        s_matrix = nt.Matrix_ps(self.overlap)
        nt.DensityMatrixSolvers.McWeenyStep(d_matrix, s_matrix, dout_matrix)
        dout_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        # Compare
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            normval = abs(norm(result_s - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        comm.barrier()

    def test_energy_density(self):
        '''Test the routines to compute the weighted-energy density matrix.'''
        from helpers import THRESHOLD, result_file
        from scipy.sparse.linalg import norm
        from scipy.io import mmread

        # Reference Solution
        self.create_matrices()
        fmat = mmread(self.hamiltonian)
        dmat = mmread(self.density)
        edm = dmat.dot(fmat).dot(dmat)

        # NTPoly
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        density_matrix = nt.Matrix_ps(self.density)
        edm_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        nt.DensityMatrixSolvers.EnergyDensityMatrix(
            fock_matrix, density_matrix, edm_matrix)

        edm_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        # Compare
        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            normval = abs(norm(edm - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        comm.barrier()


class TestChemistry_r(TestChemistry, unittest.TestCase):
    '''Specialization for real matrices'''
    # complex test
    is_complex = False

    def testrealio(self):
        '''Test routines to read data produced by a real chemistry program.'''
        from helpers import THRESHOLD, result_file
        from scipy.sparse.linalg import norm
        from scipy.io import mmread

        density_matrix = nt.Matrix_ps(self.realio)
        density_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        normval = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            self.CheckMat = mmread(self.realio)
            normval = abs(norm(self.CheckMat - ResultMat))
        global_norm = comm.bcast(normval, root=0)
        self.assertLessEqual(global_norm, THRESHOLD)
        comm.barrier()

    def test_PExtrapolate(self):
        '''Test the density extrapolation routine.'''
        from helpers import result_file

        f1 = nt.Matrix_ps(self.geomh1)
        o1 = nt.Matrix_ps(self.geomo1)
        o2 = nt.Matrix_ps(self.geomo2)
        isqm1 = nt.Matrix_ps(f1.GetActualDimension())
        d1 = nt.Matrix_ps(f1.GetActualDimension())
        extrapd = nt.Matrix_ps(f1.GetActualDimension())

        permutation = nt.Permutation(f1.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(
            o1, isqm1, self.solver_parameters)
        nt.DensityMatrixSolvers.TRS2(
            f1, isqm1, self.nel, d1, self.solver_parameters)

        nt.GeometryOptimization.PurificationExtrapolate(d1, o2, self.nel,
                                                        extrapd,
                                                        self.solver_parameters)
        extrapd.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full_extrap()
        comm.barrier()

    def test_SExtrapolate(self):
        from helpers import result_file
        '''Test the density extrapolation routine.'''
        f1 = nt.Matrix_ps(self.geomh1)
        o1 = nt.Matrix_ps(self.geomo1)
        o2 = nt.Matrix_ps(self.geomo2)
        isqm1 = nt.Matrix_ps(f1.GetActualDimension())
        d1 = nt.Matrix_ps(f1.GetActualDimension())
        extrapd = nt.Matrix_ps(f1.GetActualDimension())

        permutation = nt.Permutation(f1.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.solver_parameters.SetLoadBalance(permutation)

        nt.SquareRootSolvers.InverseSquareRoot(
            o1, isqm1, self.solver_parameters)
        nt.DensityMatrixSolvers.TRS2(
            f1, isqm1, self.nel, d1, self.solver_parameters)

        nt.GeometryOptimization.LowdinExtrapolate(d1, o1, o2, extrapd,
                                                  self.solver_parameters)
        extrapd.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_full_extrap()
        comm.barrier()

    def test_cholesky_basis(self):
        """
        Test that our solvers work even with a Cholesky basis.
        """
        from helpers import result_file
        # Setup
        self.create_matrices()
        fock_matrix = nt.Matrix_ps(self.hamiltonian)
        overlap_matrix = nt.Matrix_ps(self.overlap)
        inverse_sqrt_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        density_matrix = nt.Matrix_ps(fock_matrix.GetActualDimension())
        lmatrix = nt.Matrix_ps(fock_matrix.GetActualDimension())

        # Invert the Overlap Matrix
        nt.LinearSolvers.CholeskyDecomposition(overlap_matrix, lmatrix,
                                               self.solver_parameters)
        nt.InverseSolvers.Invert(lmatrix, inverse_sqrt_matrix,
                                 self.solver_parameters)

        for routine in [nt.DensityMatrixSolvers.PM,
                        nt.DensityMatrixSolvers.TRS2,
                        nt.DensityMatrixSolvers.TRS4,
                        nt.DensityMatrixSolvers.HPCP,
                        nt.DensityMatrixSolvers.ScaleAndFold]:
            print("Solver: ", routine)
            if routine == nt.DensityMatrixSolvers.ScaleAndFold:
                cp, homo, lumo = self.compute_cp()
                routine(fock_matrix, inverse_sqrt_matrix, self.nel,
                        density_matrix, homo, lumo, self.solver_parameters)
            else:
                routine(fock_matrix, inverse_sqrt_matrix, self.nel,
                        density_matrix, self.solver_parameters)
            density_matrix.WriteToMatrixMarket(result_file)
            comm.barrier()
            self.check_full()
            comm.barrier()


class TestChemistry_c(TestChemistry, unittest.TestCase):
    '''Specialization for complex matrices.'''
    is_complex = True


if __name__ == '__main__':
    unittest.main()
