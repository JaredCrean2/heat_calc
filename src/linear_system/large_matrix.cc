#include "linear_system/large_matrix.h"

namespace linear_system {

void solve(LargeMatrixPtr mat, const DiscVectorPtr b, DiscVectorPtr x)
{
  auto& b_vec = b->getVector();
  auto& x_vec = x->getVector();

  mat->solve(b_vec, x_vec);
}

} // namespace