'''generate.py a program to generate a random graph's exponential.

Usage:
    python generate.py number_of_nodes matrix_file exponential_file
'''
from sys import argv
from networkx import erdos_renyi_graph, to_scipy_sparse_matrix
from scipy.linalg import funm
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from numpy import diag, exp

if __name__ == "__main__":
    # Generate a random directed graph
    nodes = int(argv[1])
    prob = 0.05
    graph = erdos_renyi_graph(nodes, prob, directed=True)
    matrix = to_scipy_sparse_matrix(graph)

    # Compute The Exponential
    emat = funm(matrix.todense(), lambda x: exp(x))

    # Write To File
    mmwrite(argv[2], csr_matrix(matrix * 1.0))
    mmwrite(argv[3], csr_matrix(emat))
