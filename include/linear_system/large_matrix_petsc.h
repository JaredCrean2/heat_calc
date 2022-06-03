#ifndef LARGE_MATRIX_PETSC_H
#define LARGE_MATRIX_PETSC_H

#include "large_matrix.h"
#include "sparsity_pattern.h"
#include "petscksp.h"
#include "petscmat.h"
#include <map>

namespace linear_system {

struct LargeMatrixOptsPetsc : public LargeMatrixOpts
{
  LargeMatrixOptsPetsc()
  {
    static int i = 0;
    opts_prefix = std::string("prefix") + std::to_string(i++) + "_";
  }

  std::string opts_prefix;  // prefix to add to options when putting into Petsc database
  std::map<std::string, std::string> petsc_opts;  // options themselves (keys should *not* include opts_preifx)
};

void setPetscOptions(const LargeMatrixOptsPetsc& opts);


class LargeMatrixPetsc : public LargeMatrix
{
  public:
    LargeMatrixPetsc(DofInt mlocal, DofInt nlocal, LargeMatrixOptsPetsc opts, std::shared_ptr<SparsityPattern> sparsity_pattern);

    ~LargeMatrixPetsc();

  protected:

    void zeroMatrix_impl() override
    {
      MatZeroEntries(m_A);
      m_new_matrix_values = true;
    }

    // dofs are *global* dofs
    void assembleValues_impl(const std::vector<DofInt>& dofs, const ArrayType<Real, 2>& jac) override
    {
      MatSetValues(m_A, dofs.size(), dofs.data(), dofs.size(), dofs.data(), jac.data(), ADD_VALUES);
    }

    void finishMatrixAssembly_impl() override;

    // factor the matrix/update the preconditioner.
    void factor_impl() override;

    // vectors are *local* vectors (owned + ghost)
    // solve Ax = b for x.  The initial value of x might be used for an initial guess.
    // factor() must have been called first
    void solve_impl(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x) override;

    void matVec_impl(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override;

  private:

    void copyVec(const ArrayType<Real, 1>& vec_in, Vec vec_out, bool set_ghosts);

    void copyVec(const Vec vec_in, ArrayType<Real, 1>& vec_out, bool set_ghosts);

    LargeMatrixOptsPetsc m_opts;
    Vec m_x;
    Vec m_b;
    Mat m_A;
    KSP m_ksp;
    PC m_pc;
    bool m_new_matrix_values = true;
    std::vector<PetscInt> m_owned_dof_to_local;
    std::vector<PetscInt> m_ghost_dofs_to_local;
};


} // namespace

#endif