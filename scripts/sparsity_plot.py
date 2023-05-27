#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('AGG')
import numpy as np
import scipy
import sys
import ilupp

np.set_printoptions(threshold=sys.maxsize)


def createMatrix():
  n = 1001
  A = np.zeros((n, n))

  for i in range(n):
    if (i-1 >= 0):
      A[i, i-1] = 1

    A[i, i] = 2

    if (i+1 < n):
      A[i, i+1] = 1

    A[-1, (int(n/4)):n] = 3
    A[(int(n/4)):n, -1] = 4

  return A

def computePreconditionedOperator(A):
  Ailu = scipy.sparse.linalg.spilu(A)

  Ailu_dense = np.zeros(A.shape)
  n = A.shape[0]
  b = np.zeros(n)
  for i in range(n):
    b.fill(0)
    b[i] = 1
    Ailu_dense[i, :] = Ailu.solve(b)

  print("Ailu_dense = ", Ailu_dense)

  print("preconditioned operator = ", np.matmul(Ailu_dense, A))

  return np.matmul(Ailu_dense, A)

def computeILU0PreconditionedOperator(A):
  Asparse = scipy.sparse.csr_matrix(A)
  print("Asparse nnz = ", Asparse.getnnz())
  Ailu = ilupp.ILU0Preconditioner(Asparse)

  Ailu_dense = np.zeros(A.shape)
  n = A.shape[0]
  b = np.zeros(n)
  for i in range(n):
    b.fill(0)
    b[i] = 1
    Ailu_dense[i, :] = Ailu @ b

  print("Ail0u_dense = ", Ailu_dense)

  return np.matmul(Ailu_dense, A)



def plotEigenvalues(A, ax):
  eigs = np.sort(np.linalg.eigvals(A))

  eigs_real = np.real(eigs)
  eigs_complex = np.imag(eigs)

  ax.scatter(eigs_real, eigs_complex)
  ax.set_xlabel("Real")
  ax.set_ylabel("Imag")



A = createMatrix()

fig, axs = plt.subplots(3)

print("Ainv = ", np.linalg.inv(A))
plotEigenvalues(A, axs[0])
plotEigenvalues(np.linalg.inv(A), axs[1])
plotEigenvalues(computeILU0PreconditionedOperator(A), axs[2])
fig.savefig("eigenvalues.png", dpi=600)

