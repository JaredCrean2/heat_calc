#include "test_helper.h"

ArrayType<Real, 2> make_mat(int m, int n, const std::vector<Real>& vals)
{
  assert(vals.size() == static_cast<unsigned int>(m*n));
   ArrayType<Real, 2> arr(boost::extents[m][n]);
   for (int i=0; i < m; ++i)
     for (int j=0; j < n; ++j)
       arr[i][j] = vals[i * n + j];

  return arr;
}

ArrayType<Real, 1> make_vec(const std::vector<Real>& vals)
{
  ArrayType<Real, 1> arr(boost::extents[vals.size()]);
  for (unsigned int i=0; i < vals.size(); ++i)
   arr[i] = vals[i];

  return arr;
}

