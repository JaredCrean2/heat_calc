#include "linear_system/large_matrix.h"

namespace linear_system {

void solve(LargeMatrixPtr mat, const DiscVectorPtr b, DiscVectorPtr x)
{
  auto& b_vec = b->getVector();
  auto& x_vec = x->getVector();
  assert(b_vec.shape()[0] == mat->getMLocal());
  assert(x_vec.shape()[0] == mat->getMLocal());

  mat->solve(b_vec, x_vec);
}

} // namespace