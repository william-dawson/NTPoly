"""
Visualize two matrices side by side.
"""
import matplotlib.pyplot as plt
from scipy.io import mmread

###############################################################################
if __name__ == "__main__":
    input_mat = mmread("input.mtx")
    output_mat = mmread("output.mtx")
    plt.subplot(211)
    plt.imshow(input_mat.todense())
    plt.subplot(212)
    plt.imshow(output_mat.todense())
    plt.show()
