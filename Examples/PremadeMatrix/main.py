# Generic modules
import sys

# NTPoly
import NTPolySwig as nt

# MPI Module
from mpi4py import MPI
comm = MPI.COMM_WORLD

##########################################################################
if __name__ == "__main__":
    rank = comm.Get_rank()

    # Process The Input Parameters
    for i in range(1, len(sys.argv), 2):
        argument = sys.argv[i]
        argument_value = sys.argv[i + 1]
        if argument == '--hamiltonian':
            hamiltonian_file = argument_value
        elif argument == '--overlap':
            overlap_file = argument_value
        elif argument == '--density':
            density_file_out = argument_value
        elif argument == '--process_rows':
            process_rows = int(argument_value)
        elif argument == '--process_columns':
            process_columns = int(argument_value)
        elif argument == '--process_slices':
            process_slices = int(argument_value)
        elif argument == '--number_of_electrons':
            number_of_electrons = int(argument_value)
        elif argument == '--threshold':
            threshold = float(argument_value)
        elif argument == '--convergence_threshold':
            convergence_threshold = float(argument_value)

    # Setup the process grid.
    nt.ConstructProcessGrid(process_rows, process_columns, process_slices)

    # Read in the matrices from file.
    Hamiltonian = nt.DistributedSparseMatrix(hamiltonian_file)
    Overlap = nt.DistributedSparseMatrix(overlap_file)
    ISQOverlap = nt.DistributedSparseMatrix(Hamiltonian.GetActualDimension())
    Density = nt.DistributedSparseMatrix(Hamiltonian.GetActualDimension())
    chemical_potential = 0

    # Set Up The Solver Parameters.
    permutation = nt.Permutation(Hamiltonian.GetLogicalDimension())
    permutation.SetRandomPermutation()
    solver_parameters = nt.IterativeSolverParameters()
    solver_parameters.SetConvergeDiff(convergence_threshold)
    solver_parameters.SetThreshold(threshold)
    solver_parameters.SetLoadBalance(permutation)
    solver_parameters.SetVerbosity(True)

    # Call the solver routines.
    nt.SquareRootSolvers.InverseSquareRoot(
        Overlap, ISQOverlap, solver_parameters)
    chemical_potential = nt.DensityMatrixSolvers.TRS2(Hamiltonian, ISQOverlap,
                                 number_of_electrons,
                                 Density, solver_parameters)

    # Print the density matrix to file.
    Density.WriteToMatrixMarket(density_file_out)
