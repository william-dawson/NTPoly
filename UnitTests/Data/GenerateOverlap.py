import scipy
import scipy.io
import numpy.random
import sys
import numpy.linalg
import time
import itertools
import random
import scipy.sparse.linalg

# Example: python GenerateOverlap.py 64 0.01 22 Overlap64.mtx
if __name__ == "__main__":
  if len(sys.argv) < 5:
    raise Exception, 'Matrix dimension, sparsity, cutoff, outfile'
  dimension = int(sys.argv[1])
  sparsity = float(sys.argv[2])
  cutoff = int(sys.argv[3])
  ofile1 = sys.argv[4]
  print "Generating..."

  # Sparsity Parameters
  banded_percentage = float(cutoff*dimension)/float(dimension**2)
  sparsity = sparsity/banded_percentage
  if (sparsity > 1.0):
    raise Exception, "Impossible to fill to required sparsity"

  # Create The Matrix
  matrix = numpy.zeros((dimension,cutoff+1))
  for j in range(0,dimension):
    matrix[j,0] = 1.0
    for i in range(1,cutoff+1):
      if (i + j < dimension):
        matrix[j,i] = random.random()/((i+1)**2)

  print "Shuffling..."
  # Shuffle The Matrix
  shuffled_ordering = range(0,dimension)
  random.shuffle(shuffled_ordering)

  #shuffled_matrix = numpy.zeros((dimension,dimension))
  rows_shuffled = []
  columns_shuffled = []
  data_shuffled = []
  for j in range(0,dimension):
    rows_shuffled.append(shuffled_ordering[j])
    columns_shuffled.append(shuffled_ordering[j])
    data_shuffled.append(matrix[j,0])
    for i in range(1,cutoff+1):
      if (i + j < dimension):
        row = j
        column = i+j
        rows_shuffled.append(shuffled_ordering[row])
        columns_shuffled.append(shuffled_ordering[column])
        data_shuffled.append(matrix[j,i])
        # Symmetry
        rows_shuffled.append(shuffled_ordering[column])
        columns_shuffled.append(shuffled_ordering[row])
        data_shuffled.append(matrix[j,i])

  shuffled_matrix = scipy.sparse.csc_matrix((data_shuffled,(rows_shuffled,columns_shuffled)),shape=(dimension,dimension))
  print "Size: ", shuffled_matrix.nnz, shuffled_matrix.nnz/float(dimension**2)

  print "Inverting"
  #inv_mat = scipy.sparse.linalg.inv(shuffled_matrix)
  #product = inv_mat.dot(shuffled_matrix)
  #product >= 1e-5
  #print product
  #inv_mat >= 1e-3
  #print "Inv Size: ", inv_mat.nnz, inv_mat.nnz/float(dimension**2)

  print "Writing"
  # Write To File
  scipy.io.mmwrite(ofile1,shuffled_matrix,symmetry="general")
  #scipy.io.mmwrite(ofile2,scipy.sparse.csr_matrix(matrix2),symmetry="general")
  
