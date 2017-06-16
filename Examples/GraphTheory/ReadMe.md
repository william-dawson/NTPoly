# Graph Theory Example

In this example, we will show you how to use the library outside of a quantum
chemistry context. This is a very simple example of a graph theory application.
In this situation, we have a network of computing nodes. The nodes are
lined up next to each other, and a network wire is connected between all
nearest neighbors.

Next, we imagine trying to improve the network performance by running some
extra wires between various machines. This situation is modeled by constructing
an adjacency matrix representation of the graph. If Nodes I and J are
connected, then A_ij is nonzero. In this example, we pick some arbitrary
values for the matrix elements. Once the matrix has been set up,
we compute the resolvent of the matrix (cI - A), which is a necessary building
block for computing the Katz centrality measure.

## Build System

See the premade matrix example for details. Build with something like:
mpif90 main.f90 -o example \
  -I../../Build/include \
  -L../../Build/lib -lNTPoly -fopenmp

And run with:
mpirun -np 4 ./example \
--process_rows 2 --process_columns 2 --process_slices 1 \
--threshold 1e-6 --number_of_nodes 2048 --extra_connections 128 \
--attenuation 0.7 --output_file Output.mtx

And the python version:
export PYTHONPATH=../../Build/python
mpirun -np 4 python main.py \
--process_rows 2 --process_columns 2 --process_slices 1 \
--threshold 1e-6 --number_of_nodes 2048 --extra_connections 128 \
--attenuation 0.7 --output_file Output.mtx
