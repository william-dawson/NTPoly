"""
Compare two exponentials to see if there is a similarity in their centrality
measures.

Usage:
    python compare.py mat1 mat2
"""
from matplotlib import pyplot as plt
from numpy import diag, array
from scipy.io import mmread
from sys import argv

if __name__ == "__main__":
    mat1 = mmread(argv[1])
    mat2 = mmread(argv[2])
    diag1 = array(diag(mat1.todense()))
    diag2 = array(diag(mat2.todense()))

    fig, ax = plt.subplots(1, 1)
    ax.plot(diag1, diag2, '.')
    ax.set_xlabel("Exact")
    ax.set_ylabel("Guo Approximation")
    plt.show()
