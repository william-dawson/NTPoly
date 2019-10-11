"""
Visualize a matrix with a log scale.
"""
from sys import argv
from scipy.io import mmread
from numpy import seterr, floor, abs, log10
from matplotlib import pyplot as plt


if __name__ == "__main__":
    # This will allow us to catch floating point errors ourself.
    seterr(all='raise')

    # Get the input matrix
    file_1 = argv[1]
    matrix1 = mmread(file_1)

    # Convert to dense so we can set everything to a log scale
    printable = matrix1.todense()
    for j in range(0, printable.shape[1]):
        for i in range(0, printable.shape[0]):
            try:
                printable[j, i] = floor(log10(abs(printable[j, i])))
            except FloatingPointError:
                printable[j, i] = -10

    # Plot
    fig, ax = plt.subplots()
    cax = ax.matshow(printable)
    cbar = fig.colorbar(cax)

    plt.show()
