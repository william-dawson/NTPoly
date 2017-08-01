################################################################################
import matplotlib.pyplot as plt
import numpy
from matplotlib.colors import LogNorm

# A function that reads in matrix market files produced by the library.
# @param file_name file to read.
def ReadMM(file_name):
    matrix = numpy.zeros((0, 0))
    with open(file_name, "r") as fin:
        line = ""
        while True:
            line = fin.next()
            if line.lstrip()[0] == "%":
                continue
            else:
                break
        split = line.split()
        matrix = numpy.zeros((int(split[0]), int(split[1])))
        for line in fin:
            split = line.split()
            if (len(split) == 3):
                matrix[int(split[0]) - 1, int(split[1]) - 1] = float(split[2])
    return matrix


################################################################################
if __name__ == "__main__":
    input_mat = ReadMM("input.mtx")
    output_mat = ReadMM("output.mtx")
    plt.subplot(211)
    plt.imshow(input_mat)
    plt.subplot(212)
    plt.imshow(output_mat)
    plt.show()
