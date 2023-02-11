#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include "mpi.h"

int commRank(MPI_Comm comm);

int commSize(MPI_Comm comm);

#endif