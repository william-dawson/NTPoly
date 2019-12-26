"""
Visualize a matrix.
"""
import matplotlib.pyplot as plt
from scipy.io import mmread

###############################################################################
if __name__ == "__main__":
    density_mat = mmread("Density.mtx")
    plt.imshow(density_mat.todense())
    plt.show()
