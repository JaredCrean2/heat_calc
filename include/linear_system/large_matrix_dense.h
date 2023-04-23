#ifndef LARGE_MATRIX_DENSE_H
#define LARGE_MATRIX_DENSE_H

#include "large_matrix.h"
#include "bla_wrapper.h"
#include "linear_system/sparsity_pattern_dense.h"

#include <ostream>

namespace linear_system {

class LargeMatrixDense : public LargeMatrix
{
  public:
    LargeMatrixDense(DofInt mlocal, DofInt nlocal, LargeMatrixOpts opts) :
      LargeMatrix(mlocal, nlocal, std::make_shared<SparsityPatternDense>(mlocal)),
      m_opts(opts),
      m_matrix(mlocal * nlocal)
    {}

    const Real& operator()(int i, int j) const { return m_matrix[getIdx(i, j)]; }

    void printToStdout();

  protected:
    void zeroMatrix_impl() override;

    void assembleValues_impl(const std::vector<DofInt>& dofs_rows, const std::vector<DofInt>& dofs_cols, const ArrayType<Real, 2>& jac) override;

    void factor_impl() override;

    void solve_impl(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x) override;

    void matVec_impl(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override;


  private:
    DofInt getIdx(DofInt m, DofInt n) const
    {
      assert(m >= 0 && m < getMLocal());
      assert(n >= 0 && n < getNLocal());
 
      // note: this is column major to make Lapack interface easier
      return m + n * getMLocal();
    }

    LargeMatrixOpts m_opts;
    std::vector<Real> m_matrix;
    std::vector<Real> m_matrix_factorization;
    std::vector<lapack_int> m_ipiv;

    friend std::ostream& operator<<(std::ostream& os, const LargeMatrixDense& mat);

};

std::ostream& operator<<(std::ostream& os, const LargeMatrixDense& mat);


} // namespace

#endif