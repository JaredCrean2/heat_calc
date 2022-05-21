#ifndef TEST_HELPER_H
#define TEST_HELPER_H

#include "ProjectDefs.h"
#include <vector>
#include "mpi.h"

ArrayType<Real, 2> make_mat(int m, int n, const std::vector<Real>& vals);

ArrayType<Real, 1> make_vec(const std::vector<Real>& vals);



#define SERIAL_ONLY()                              \
{                                                  \
  int comm_size;                                   \
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);       \
  if (comm_size > 1)                               \
    GTEST_SKIP();                                  \
}                                                  \

#define SERIAL_ONLY_RETURN()                              \
{                                                  \
  int comm_size;                                   \
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);       \
  if (comm_size > 1)                               \
    return;                                        \
} 


#endif
