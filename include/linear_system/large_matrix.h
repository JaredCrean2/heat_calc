#ifndef LARGE_MATRIX_H
#define LARGE_MATRIX_H

#include "ProjectDefs.h"
#include "discretization/disc_vector.h"
#include "linear_system/sparsity_pattern.h"
#include "utils/error_handling.h"

namespace linear_system {

struct LargeMatrixOpts
{
  virtual ~LargeMatrixOpts() = default;

  bool is_structurally_symmetric = false;
  bool is_value_symmetric        = false;
  bool factor_in_place           = false;
};


class LargeMatrix
{
  public:
    LargeMatrix(DofInt m_local, DofInt n_local, std::shared_ptr<SparsityPattern> sparsity) :
      m_mlocal(m_local),
      m_nlocal(n_local),
      m_sparsity_pattern(sparsity),
      m_is_factored(false)
    {}

    virtual ~LargeMatrix() = default;

    LargeMatrix(const LargeMatrix&) = delete;

    LargeMatrix& operator=(const LargeMatrix&) = delete;

    DofInt getMLocal() const { return m_mlocal; }

    DofInt getNLocal() const { return m_nlocal; }

    std::array<DofInt, 2> getSize() const { return {m_mlocal, m_nlocal}; }

    std::shared_ptr<SparsityPattern> getSparsityPattern() const { return m_sparsity_pattern; }

    void zeroMatrix()
    {
      m_is_factored = false;
      zeroMatrix_impl();
    }

    // if dof < 0, the corresponding entries are ignored
    void assembleValues(const std::vector<DofInt>& dofs, const ArrayType<Real, 2>& jac)
    {
#ifndef NDEBUG
      assert(!m_is_factored);

      assert(jac.shape()[0] == dofs.size());
      assert(jac.shape()[1] == dofs.size());

      // check uniqueness
      std::vector<DofInt> dofs_copy;
      for (unsigned int i=0; i < dofs.size(); ++i)
        if (dofs[i] >= 0)
          dofs_copy.push_back(i);

      auto it = std::unique(dofs_copy.begin(), dofs_copy.end());
      assert(it == dofs_copy.end());
#endif

      assembleValues(dofs, dofs, jac);
    }

    void assembleValues(const std::vector<DofInt>& dofs_rows, const std::vector<DofInt>& dofs_cols, const ArrayType<Real, 2>& jac)
    {
      assert(!m_is_factored);

      assert(jac.shape()[0] == dofs_rows.size());
      assert(jac.shape()[1] == dofs_cols.size());    

      assembleValues_impl(dofs_rows, dofs_cols, jac);  

    }
        
    void finishMatrixAssembly() { finishMatrixAssembly_impl(); };


    // factor the matrix/update the preconditioner.
    void factor()
    {
      factor_impl();
      m_is_factored = true;
    }

    // solve Ax = b for x.  The initial value of x might be used for an initial guess.
    // factor() must have been called first
    void solve(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x)
    {
      assertAlways(m_is_factored, "must factor matrix before solving a linear system");
      solve_impl(b, x);
    }

    // compute b = A * x, overwriting b
    void matVec(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
    {
      assert(x.shape()[0] >= getNLocal());
      assert(b.shape()[0] >= getMLocal());
      matVec_impl(x, b);
    }

  protected:

    bool getIsFactored() const { return m_is_factored; }

    virtual void zeroMatrix_impl() = 0;

    virtual void assembleValues_impl(const std::vector<DofInt>& dofs_rows, const std::vector<DofInt>& dofs_cols, const ArrayType<Real, 2>& jac) = 0;

    virtual void finishMatrixAssembly_impl() {};

    // factor the matrix/update the preconditioner.
    virtual void factor_impl() = 0;

    // solve Ax = b for x.  The initial value of x might be used for an initial guess.
    // factor() must have been called first
    virtual void solve_impl(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x) = 0;

    virtual void matVec_impl(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) = 0;

  private:
    DofInt m_mlocal;
    DofInt m_nlocal;
    std::shared_ptr<SparsityPattern> m_sparsity_pattern;
    bool m_is_factored;
};

using LargeMatrixPtr = std::shared_ptr<LargeMatrix>;

void solve(LargeMatrixPtr mat, const DiscVectorPtr b, DiscVectorPtr x);

} // namespace

#endif
