#ifndef LARGE_MATRIX_PETSC_H
#define LARGE_MATRIX_PETSC_H

#include "large_matrix.h"
#include "sparsity_pattern.h"
#include "petscksp.h"
#include "petscmat.h"
#include <map>
#include <petscviewer.h>

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

// sets an option in the Petsc global options database.  key should *not* have the
// hyphen prefix.
void setPetscGlobalOption(const std::string& key, const std::string& val);


class LargeMatrixPetsc : public LargeMatrix
{
  public:
    LargeMatrixPetsc(LargeMatrixOptsPetsc opts, std::shared_ptr<SparsityPattern> sparsity_pattern);

    ~LargeMatrixPetsc();

    std::shared_ptr<LargeMatrix> clone() override;

    void printToStdout()
    {
      MatView(m_A, PETSC_VIEWER_STDOUT_WORLD);
    }  

    LargeMatrixPetsc(std::shared_ptr<SparsityPattern> sparsity_pattern);

  protected:

    void zeroMatrix_impl() override
    {
      MatZeroEntries(m_A);
      m_new_matrix_values = true;
    }

    // dofs are *global* dofs
    void assembleValues_impl(const std::vector<DofInt>& dofs_rows, const std::vector<DofInt>& dofs_cols, const ArrayType<Real, 2>& jac) override
    {
      MatSetValues(m_A, dofs_rows.size(), dofs_rows.data(), dofs_cols.size(), dofs_cols.data(), jac.data(), ADD_VALUES);
    }

    void finishMatrixAssembly_impl() override;

    // factor the matrix/update the preconditioner.
    void factor_impl() override;

    // vectors are *local* vectors (owned + ghost)
    // solve Ax = b for x.  The initial value of x might be used for an initial guess.
    // factor() must have been called first
    void solve_impl(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x) override;

    void matVec_impl(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override;

    void axpy_impl(Real alpha, std::shared_ptr<LargeMatrix> x) override;  

  private:

    void copyVec(const ArrayType<Real, 1>& vec_in, Vec vec_out, bool set_ghosts);

    void copyVec(const Vec vec_in, ArrayType<Real, 1>& vec_out, bool set_ghosts);

    LargeMatrixOptsPetsc m_opts;
    Vec m_x = nullptr;
    Vec m_b = nullptr;
    Mat m_A = nullptr;
    KSP m_ksp = nullptr;
    PC m_pc = nullptr;
    bool m_new_matrix_values = true;
    std::vector<PetscInt> m_owned_dof_to_local;
    std::vector<PetscInt> m_ghost_dofs_to_local;
};


} // namespace

#endif