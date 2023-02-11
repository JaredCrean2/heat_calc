#include "utils/mpi_utils.h"

int commRank(MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  return rank;
}

int commSize(MPI_Comm comm)
{
  int size;
  MPI_Comm_size(comm, &size);

  return size;
}