# NTPoly
import NTPolySwig as nt

# MPI Module
from mpi4py import MPI
comm = MPI.COMM_WORLD


def ConstructGuoMatrix(InMat, OutMat):
    '''
    Construct the Hermitian matrix from a nonsymmetric matrix.
    '''
    # First Symmetrize The Input Matrix.
    tlist = nt.TripletList_r()
    stlist = nt.TripletList_r()
    InMat.GetTripletList(tlist)

    for i in range(0, tlist.GetSize()):
        temp = tlist.GetTripletAt(i)
        stlist.Append(temp)
        if temp.index_row != temp.index_column:
            temp2 = nt.Triplet_r()
            temp2.index_row = temp.index_column
            temp2.index_column = temp.index_row
            temp2.point_value = temp.point_value
            stlist.Append(temp2)
    SMat = nt.Matrix_ps(InMat.GetActualDimension())
    SMat.FillFromTripletList(stlist)

    # Construct The Guide Matrix.
    Guide = nt.Matrix_ps(InMat.GetActualDimension())
    Guide.Increment(InMat, -1.0)

    # Now iterate over the entries in the guide matrix.
    clist = nt.TripletList_c()
    Guide.GetTripletList(clist)
    temp_c = nt.Triplet_c()
    for i in range(0, clist.GetSize()):
        temp = tlist.GetTripletAt(i)
        temp_c.index_row = temp.index_row
        temp_c.index_column = temp.index_column
        temp_c.point_value = 0.0 + 1j
        clist.Append(temp_c)
    CMatrix = nt.Matrix_ps(InMat.GetActualDimension())
    CMatrix.FillFromTripletList(clist)

    OutMat.Transpose(CMatrix)
    OutMat.Conjugate()
    OutMat.Increment(CMatrix)
    OutMat.Increment(SMat)


##########################################################################
if __name__ == "__main__":
    from sys import argv

    rank = comm.Get_rank()
    total_processors = comm.Get_size()

    # Process The Input Parameters
    for i in range(1, len(argv), 2):
        argument = argv[i]
        argument_value = argv[i + 1]
        if argument == '--input_file':
            input_file = argument_value
        if argument == '--exponential_file':
            exponential_file = argument_value
        elif argument == '--process_rows':
            process_rows = int(argument_value)
        elif argument == '--process_columns':
            process_columns = int(argument_value)
        elif argument == '--process_slices':
            process_slices = int(argument_value)
        elif argument == '--threshold':
            threshold = float(argument_value)

    # Setup the process grid.
    nt.ConstructGlobalProcessGrid(
        process_rows, process_columns, process_slices)

    # Construct The Hermitian Matrix
    InMat = nt.Matrix_ps(input_file)
    GMat = nt.Matrix_ps(InMat.GetActualDimension())
    ConstructGuoMatrix(InMat, GMat)
    GMat.Scale(0.5)

    # Compute The Exponential
    solver_parameters = nt.SolverParameters()
    solver_parameters.SetThreshold(threshold)
    OMat = nt.Matrix_ps(InMat.GetActualDimension())
    nt.ExponentialSolvers.ComputeExponential(GMat, OMat, solver_parameters)

    # Write To File
    OMat.WriteToMatrixMarket(exponential_file)
    nt.DestructGlobalProcessGrid()
