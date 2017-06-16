##########################################################################
import sys
import scipy
import scipy.io
import numpy
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

##########################################################################
if __name__ == "__main__":
    numpy.seterr(all='raise')
    file_1 = sys.argv[1]

    matrix1 = scipy.io.mmread(file_1)
    printable = matrix1.todense()
    for j in range(0, printable.shape[1]):
        for i in range(0, printable.shape[0]):
            try:
                printable[j, i] = numpy.floor(
                    numpy.log10(numpy.abs(printable[j, i])))
            except:
                printable[j, i] = -10

    fig, ax = plt.subplots()
    cax = ax.matshow(printable)
    cbar = fig.colorbar(cax)

    plt.show()
