#ifndef HEAT_CALC_PROJECT_DEFS_H
#define HEAT_CALC_PROJECT_DEFS_H

#include "boost/multi_array.hpp"
#include "mpi.h"

template <typename T, long unsigned N>
using ArrayType = boost::multi_array<T, N>;
using range = boost::multi_array_types::index_range;

using Real = double;
constexpr MPI_Datatype REAL_MPI_DATATYPE = MPI_DOUBLE;

constexpr Real MIN_REAL = std::numeric_limits<Real>::lowest();
constexpr Real MAX_REAL = std::numeric_limits<Real>::max();

using DofInt = int;      // for degrees of freedom
using Index = int;       // for element/face numbers

constexpr MPI_Datatype DofInt_MPI_DATATYPE = MPI_INT;
using LocalIndex = int;  // for local to the element numbering

//TODO: make a generated header
#define MESH_USE_MDS_NUMBERING

template <typename T>
void fill(ArrayType<T, 1>& u, T val)
{
  for (auto& v : u)
    v = val;
}

template <typename T, typename T2>
void fill(ArrayType<T, 1>& u, T2 val)
{
  fill(u, T(val));
}
/*
//TODO: there is another way to do this with variadic templates
template <typename Array, typename T1>
view(Array& A, const T& t1) -> decltype(A[boost::indices[t1]])
{
  return A[boost::indices[t1]];
}

template <typename Array, typename T1, typename T2>
view(Array& A, const T1& t1, const T2& t2) -> decltype(A[boost::indices[t1][t2]])
{
  return A[boost::indices[t1][t2]];
}

template <typename Array, typename T1, typename T2, typename T3>
view(Array& A, const T1& t1, const T2& t2, const T3& t3) -> decltype(A[boost::indices[t1][t2][t3]])
{
  return A[boost::indices[t1][t2][t3]];
}

template <typename Array, typename T1, typename T2, typename T3, typename T4>
view(Array& A, const T1& t1, const T2& t2, const T3& t3, const T4& t4) -> decltype(A[boost::indices[t1][t2][t3][t4]])
{
  return A[boost::indices[t1][t2][t3][t4]];
}

*/
#endif