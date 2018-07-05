# Generic modules
import sys

# NT Poly
import NTPolySwig as nt

# MPI Module
from mpi4py import MPI
comm = MPI.COMM_WORLD

import numpy

import random

##########################################################################
if __name__ == "__main__":
    rank = comm.Get_rank()
    total_processors = comm.Get_size()

    # Process The Input Parameters
    for i in range(1, len(sys.argv), 2):
        argument = sys.argv[i]
        argument_value = sys.argv[i + 1]
        if argument == '--output_file':
            output_file = argument_value
        elif argument == '--process_rows':
            process_rows = int(argument_value)
        elif argument == '--process_columns':
            process_columns = int(argument_value)
        elif argument == '--process_slices':
            process_slices = int(argument_value)
        elif argument == '--threshold':
            threshold = float(argument_value)
        elif argument == '--convergence_threshold':
            convergence_threshold = float(argument_value)
        elif argument == '--attenuation':
            attenuation = float(argument_value)
        elif argument == '--number_of_nodes':
            number_of_nodes = int(argument_value)
        elif argument == '--extra_connections':
            extra_connections = int(argument_value)

    # Setup the process grid.
    nt.ConstructProcessGrid(process_rows, process_columns, process_slices)

    # Set Up The Solver Parameters.
    solver_parameters = nt.IterativeSolverParameters()
    solver_parameters.SetThreshold(threshold)
    solver_parameters.SetConvergeDiff(convergence_threshold)
    solver_parameters.SetVerbosity(True)

    # Divide The Nodes Amongst Processors.
    number_of_local_nodes = int(number_of_nodes / total_processors)
    starting_node = number_of_local_nodes * rank
    # Handles the edge case
    if rank == total_processors - 1:
        number_of_local_nodes = number_of_nodes - rank * number_of_local_nodes
    local_nodes = []
    for counter in range(0, number_of_local_nodes):
        local_nodes.append(starting_node + (counter))
    ending_node = local_nodes[-1]

    # Fill the matrix
    NetworkMat = nt.Matrix_ps(number_of_nodes)
    ResultMat = nt.Matrix_ps(number_of_nodes)

    triplet_list = nt.TripletList_r(0)
    temp_triplet = nt.Triplet_r()

    # First add the connection between each node and itself.
    for counter in range(0, number_of_local_nodes):
        temp_triplet.index_row = local_nodes[counter] + 1
        temp_triplet.index_column = local_nodes[counter] + 1
        temp_triplet.point_value = 1
        triplet_list.Append(temp_triplet)

    # Now connections between nearest neighbors.
    for counter in range(0, number_of_local_nodes):
        temp_triplet.index_row = local_nodes[counter] + 1
        temp_triplet.point_value = 0.1
        if local_nodes[counter] == 0:
            # Right value
            temp_triplet.index_column = local_nodes[counter] + 1 + 1
            triplet_list.Append(temp_triplet)
        elif local_nodes[counter] == number_of_nodes - 1:
            # Left value
            temp_triplet.index_column = local_nodes[counter] - 1 + 1
            triplet_list.Append(temp_triplet)
        else:
            # Left value
            temp_triplet.index_column = local_nodes[counter] - 1 + 1
            triplet_list.Append(temp_triplet)
            # Right value
            temp_triplet.index_column = local_nodes[counter] + 1 + 1
            triplet_list.Append(temp_triplet)

    # Finally the random extra connections.
    extra_scratch = numpy.zeros(number_of_nodes)
    counter = 0
    while counter < extra_connections:
        extra_source_node = random.randint(0, number_of_nodes-1)
        extra_destination_node = random.randint(0, number_of_nodes-1)
        if extra_scratch[extra_source_node] != 1 and \
                extra_scratch[extra_destination_node] != 1 and \
                extra_source_node != extra_destination_node and \
                extra_source_node != extra_destination_node - 1 and \
                extra_source_node != extra_destination_node + 1:
            counter = counter + 1
            extra_scratch[extra_source_node] = 1
            extra_scratch[extra_destination_node] = 1

            if extra_source_node >= starting_node and \
                    extra_source_node <= ending_node:
                temp_triplet.index_row = extra_source_node+1
                temp_triplet.index_column = extra_destination_node+1
                temp_triplet.point_value = 0.1
                triplet_list.Append(temp_triplet)
            elif extra_destination_node >= starting_node and \
                    extra_destination_node <= ending_node:
                temp_triplet.index_row = extra_destination_node+1
                temp_triplet.index_column = extra_source_node+1
                temp_triplet.point_value = 0.1
                triplet_list.Append(temp_triplet)

    # Finally build the matrix
    NetworkMat.FillFromTripletList(triplet_list)

    # Solve
    ResMat = nt.Matrix_ps(number_of_nodes)
    ResMat.FillIdentity()
    ResMat.Increment(NetworkMat,alpha=-1.0*attenuation)

    nt.InverseSolvers.Invert(ResMat,ResultMat,solver_parameters)

    # Print the density matrix to file.
    ResultMat.WriteToMatrixMarket(output_file)
