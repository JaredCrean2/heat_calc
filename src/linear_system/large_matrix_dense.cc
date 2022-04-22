#include "linear_system/large_matrix_dense.h"
#include "linear_system/bla_wrapper.h"

namespace linear_system {

void LargeMatrixDense::zeroMatrix_impl()
{
  std::fill(m_matrix.begin(), m_matrix.end(), 0);
}

void LargeMatrixDense::assembleValues_impl(const std::vector<DofInt>& dofs, const ArrayType<Real, 2>& jac)
{
  unsigned int nvals = dofs.size();
  for (int j=0; j < nvals; ++j)
  {
    if (dofs[j] < 0)
      continue;

    for (int i=0; i < nvals; ++i)
      if (dofs[i] >= 0)
        m_matrix[getIdx(dofs[i], dofs[j])] += jac[i][j];
  }
}

void LargeMatrixDense::factor_impl()
{
  
  std::vector<Real>* matrix_factorization;
  if (m_opts.factor_in_place)
  {
    matrix_factorization = &m_matrix;
  } else
  {
    matrix_factorization = &m_matrix_factorization;
    *matrix_factorization = m_matrix;
  }

  if (m_opts.is_value_symmetric)
  {
    potrf('L', getMLocal(), *matrix_factorization);
  } else
  {
    m_ipiv.resize(getMLocal());
    getrf(getMLocal(), getNLocal(), matrix_factorization->data(), getMLocal(), m_ipiv.data());
  }

}



void LargeMatrixDense::solve_impl(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x)
{
  assert(b.shape()[0] == getMLocal());
  assert(x.shape()[0] == getMLocal());

  std::copy(&(b[0]), (&b[0] + getMLocal()), &(x[0]));

  auto& matrix_factorization = m_opts.factor_in_place ? m_matrix : m_matrix_factorization;
  if (m_opts.is_value_symmetric)
  {
    potrs('L', getMLocal(), 1, matrix_factorization.data(), getMLocal(), &(x[0]), getMLocal());
  } else
  {
    getrs('N', getMLocal(), 1, matrix_factorization.data(), getMLocal(), m_ipiv.data(), &(x[0]), getMLocal());
  }  
}


void LargeMatrixDense::matVec_impl(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  assertAlways(!getIsFactored(), "cannot do mat-vec when matrix is factored");

  if (m_opts.is_value_symmetric)
    symv('L', getMLocal(), 1, m_matrix.data(), getMLocal(), &(x[0]), 1, 0, &(b[0]), 1);
  else
    gemv('N', getMLocal(), getNLocal(), 1, m_matrix.data(), getMLocal(), &(x[0]), 1, 0, &(b[0]), 1);
}



std::ostream& operator<<(std::ostream& os, const LargeMatrixDense& mat)
{
  for (DofInt i=0; i < mat.getMLocal(); ++i)
  {
    for (DofInt j=0; j < mat.getNLocal(); ++j)
    {
      os << mat.m_matrix[mat.getIdx(i, j)];
      if (j != mat.getNLocal() - 1)
        os << ", ";
    }
    if (i != mat.getMLocal() - 1)
      os << std::endl;
  }

  return os;
}


}  // namespace