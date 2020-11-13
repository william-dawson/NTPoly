# NTPoly
import NTPolySwig as nt

# MPI Module
from mpi4py import MPI
comm = MPI.COMM_WORLD

##########################################################################
if __name__ == "__main__":
    from sys import argv

    rank = comm.Get_rank()

    # Process The Input Parameters
    for i in range(1, len(argv), 2):
        argument = argv[i]
        argument_value = argv[i + 1]
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
        elif argument == '--converge_overlap':
            converge_overlap = float(argument_value)
        elif argument == '--converge_density':
            converge_density = float(argument_value)

    # Setup the process grid.
    nt.ConstructGlobalProcessGrid(
        process_rows, process_columns, process_slices)
    if nt.GetGlobalIsRoot():
        nt.ActivateLogger()

    # Read in the matrices from file.
    Hamiltonian = nt.Matrix_ps(hamiltonian_file)
    Overlap = nt.Matrix_ps(overlap_file)
    ISQOverlap = nt.Matrix_ps(Hamiltonian.GetActualDimension())
    Density = nt.Matrix_ps(Hamiltonian.GetActualDimension())
    chemical_potential = 0

    # Set Up The Solver Parameters.
    permutation = nt.Permutation(Hamiltonian.GetLogicalDimension())
    permutation.SetRandomPermutation()
    solver_parameters = nt.SolverParameters()
    solver_parameters.SetConvergeDiff(converge_overlap)
    solver_parameters.SetThreshold(threshold)
    solver_parameters.SetLoadBalance(permutation)
    solver_parameters.SetVerbosity(True)

    # Call the solver routines.
    nt.SquareRootSolvers.InverseSquareRoot(
        Overlap, ISQOverlap, solver_parameters)

    # Change the solver variable for computing the density matrix.
    solver_parameters.SetConvergeDiff(converge_density)

    # Compute the density matrix.
    energy_value, chemical_potential = \
        nt.DensityMatrixSolvers.TRS2(Hamiltonian, ISQOverlap,
                                     number_of_electrons,
                                     Density,
                                     solver_parameters)

    # Print the density matrix to file.
    Density.WriteToMatrixMarket(density_file_out)

    # Cleanup
    if nt.GetGlobalIsRoot():
        nt.DeactivateLogger()
    nt.DestructGlobalProcessGrid()
