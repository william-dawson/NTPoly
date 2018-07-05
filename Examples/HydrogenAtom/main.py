# Generic modules
import sys

# NT Poly
import NTPolySwig as nt

# MPI Module
from mpi4py import MPI
comm = MPI.COMM_WORLD

import numpy

##########################################################################
if __name__ == "__main__":
    rank = comm.Get_rank()
    total_processors = comm.Get_size()

    x_start = -6.28
    x_end = 6.28

    # Process The Input Parameters
    for i in range(1, len(sys.argv), 2):
        argument = sys.argv[i]
        argument_value = sys.argv[i + 1]
        if argument == '--convergence_threshold':
            convergence_threshold = float(argument_value)
        elif argument == '--density':
            density_file_out = argument_value
        elif argument == '--grid_points':
            grid_points = int(argument_value)
        elif argument == '--process_rows':
            process_rows = int(argument_value)
        elif argument == '--process_columns':
            process_columns = int(argument_value)
        elif argument == '--process_slices':
            process_slices = int(argument_value)
        elif argument == '--threshold':
            threshold = float(argument_value)

    # Setup the process grid.
    nt.ConstructProcessGrid(process_rows, process_columns, process_slices)

    # Set Up The Solver Parameters.
    solver_parameters = nt.SolverParameters()
    solver_parameters.SetConvergeDiff(convergence_threshold)
    solver_parameters.SetThreshold(threshold)
    solver_parameters.SetVerbosity(True)

    # Divide The Work Amongst Processors.
    local_grid_points = int(grid_points / total_processors)
    start_row = local_grid_points * rank
    if rank == total_processors - 1:
        local_grid_points = grid_points - rank * local_grid_points
    full_range = numpy.arange(grid_points)
    local_rows = full_range[start_row:start_row + local_grid_points]

    # Construct A Linear Space.
    full_x, grid_spacing = numpy.linspace(x_start, x_end, num=grid_points,
                                          retstep=True)
    x_values = full_x[start_row:start_row + local_grid_points]

    # Construct The Kinetic Energy Operator.
    triplet_list = nt.TripletList_r()

    insert_location = 0
    temp_value = nt.Triplet_r()
    for row_value in local_rows:
        temp_value.index_row = int(row_value + 1)
        # Stencil point 1.
        if row_value > 1:
            temp_value.index_column = temp_value.index_row - 2
            temp_value.point_value = (-0.5) * (-1.0 / (12.0 * grid_spacing**2))
            triplet_list.Append(temp_value)
        # Stencil point 2.
        if row_value > 0:
            temp_value.index_column = temp_value.index_row - 1
            temp_value.point_value = (-0.5) * (16.0 / (12.0 * grid_spacing**2))
            triplet_list.Append(temp_value)
        # Stencil point 3.
        temp_value.index_column = temp_value.index_row
        temp_value.point_value = (-0.5) * (-30.0 / (12.0 * grid_spacing**2))
        triplet_list.Append(temp_value)
        # Stencil point 4.
        if row_value + 1 < grid_points:
            temp_value.index_column = temp_value.index_row + 1
            temp_value.point_value = (-0.5) * (16.0 / (12.0 * grid_spacing**2))
            triplet_list.Append(temp_value)
        # Stencil point 5.
        if row_value + 2 < grid_points:
            temp_value.index_column = temp_value.index_row + 2
            temp_value.point_value = (-0.5) * (-1.0 / (12.0 * grid_spacing**2))
            triplet_list.Append(temp_value)

    KineticEnergy = nt.Matrix_ps(grid_points)
    KineticEnergy.FillFromTripletList(triplet_list)

    # Construct The Potential Energy Operator.
    triplet_list = nt.TripletList_r()
    for row_value, grid_value in zip(range(0, local_grid_points), x_values):
        temp_value.index_row = start_row + row_value + 1
        temp_value.index_column = start_row + row_value + 1
        temp_value.point_value = -1.0 / abs(grid_value)
        triplet_list.Append(temp_value)
    PotentialEnergy = nt.Matrix_ps(grid_points)
    PotentialEnergy.FillFromTripletList(triplet_list)

    # Construct The Full Hamiltonian.
    Hamiltonian = nt.Matrix_ps(KineticEnergy)
    Hamiltonian.Increment(PotentialEnergy)

    # Overlap Matrix is just the identity.
    Identity = nt.Matrix_ps(grid_points)
    Identity.FillIdentity()

    # Call the solver routine.
    Density = nt.Matrix_ps(grid_points)
    chemical_potential = nt.DensityMatrixSolvers.TRS2(Hamiltonian, Identity, 2,
                                                      Density,
                                                      solver_parameters)

    # Print the density matrix to file.
    Density.WriteToMatrixMarket(density_file_out)
