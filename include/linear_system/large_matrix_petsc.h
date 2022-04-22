#ifndef LARGE_MATRIX_PETSC_H
#define LARGE_MATRIX_PETSC_H

#include "large_matrix.h"
#include "sparsity_pattern.h"
#include "petscksp.h"
#include <petscmat.h>

namespace linear_system {

struct LargeMatrixOptsPetsc : public LargeMatrixOpts
{
  LargeMatrixOptsPetsc()
  {
    static int i;
    opts_prefix = std::string("prefix") + std::to_string(i++);
  }

  std::string opts_prefix;  // prefix to add to options when putting into Petsc database
  std::map<std::string, std::string> petsc_opts;  // options themselves (keys should *not* include opts_preifx)
};

void setPetscOptions(const LargeMatrixOptsPetsc& opts);


class LargeMatrixPetsc : public LargeMatrix
{
  public:
    LargeMatrixPetsc(DofInt mlocal, DofInt nlocal, LargeMatrixOptsPetsc opts, std::shared_ptr<SparsityPattern> sparsity_pattern) :
      LargeMatrix(mlocal, nlocal),
      m_opts(opts)
    {
      setPetscOptions(opts);

      VecCreate(PETSC_COMM_WORLD, &m_x);
      PetscObjectSetName((PetscObject)m_x, (opts.opts_prefix + "_Solution").c_str());
      VecSetOptionsPrefix(m_x, opts.opts_prefix.c_str());
      VecSetSizes(m_x, getNLocal(), PETSC_DECIDE);
      VecSetFromOptions(m_x);

      VecCreate(PETSC_COMM_WORLD, &m_b);
      PetscObjectSetName((PetscObject)m_b, (opts.opts_prefix + "_rhs").c_str());
      VecSetOptionsPrefix(m_b, opts.opts_prefix.c_str());
      VecSetSizes(m_b, getMLocal(), PETSC_DECIDE);
      VecSetFromOptions(m_b);

      MatCreate(PETSC_COMM_WORLD, &m_A);
      PetscObjectSetName((PetscObject)m_A, (opts.opts_prefix + "_rhs").c_str());
      MatSetSizes(m_A, getMLocal(), getNLocal(), PETSC_DECIDE, PETSC_DECIDE);
      MatSetOptionsPrefix(m_A, opts.opts_prefix.c_str());
      MatSetFromOptions(m_A);
      if (opts.is_value_symmetric)
        MatSetOption(m_A, MAT_SYMMETRIC, PETSC_TRUE);
      else if (opts.is_structurally_symmetric)
        MatSetOption(m_A, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);

      if (opts.is_value_symmetric || opts.is_structurally_symmetric)
        MatSetOption(m_A, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);

      MatSetOption(m_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(m_A, MAT_NEW_NONZERO_LOCATIONS,      PETSC_FALSE);
      MatSetOption(m_A, MAT_NEW_NONZERO_LOCATION_ERR,   PETSC_TRUE);
      // TODO: MAT_SORTED_FULL for preallocation
      // TODO: MatSetOption(m_A, MAT_ROWS_SORTED) and similarly for MAT_COLUMNS_SORTED (see page 56 of manual)

      MatXAIJSetPreallocation(m_A, 1, sparsity_pattern->getDiagonalCounts().data(), sparsity_pattern->getOffProcCounts().data(),
                                      sparsity_pattern->getDiagonalCountsSym().data(), sparsity_pattern->getOffProcCountsSym().data());
      

      KSPCreate(PETSC_COMM_WORLD, &m_ksp);
      PetscObjectSetName((PetscObject)m_ksp, (opts.opts_prefix + "_ksp").c_str());
      KSPSetOperators(m_ksp, m_A, m_A);
      KSPSetOptionsPrefix(m_ksp, opts.opts_prefix.c_str());
      KSPSetFromOptions(m_ksp);

      KSPGetPC(m_ksp, &m_pc);
      PetscObjectSetName((PetscObject)m_pc, (opts.opts_prefix + "_pc").c_str());
      PCSetOptionsPrefix(m_pc, opts.opts_prefix.c_str());
      PCSetFromOptions(m_pc);

      if (opts.factor_in_place)
        PCFactorSetUseInPlace(m_pc, PETSC_TRUE);
    }

    ~LargeMatrixPetsc()
    {
      VecDestroy(&m_x);
      VecDestroy(&m_b);
      MatDestroy(&m_A);
      KSPDestroy(&m_ksp);
      // KSP manages the PC, no need to destroy it explicitly
    }

  protected:

    void zeroMatrix_impl() override
    {
      MatZeroEntries(m_A);
      m_new_matrix_values = true;
    }

    void assembleValues_impl(const std::vector<DofInt>& dofs, const ArrayType<Real, 2>& jac) override
    {

      MatSetValues(m_A, dofs.size(), dofs.data(), dofs.size(), dofs.data(), jac.data(), ADD_VALUES);
    }
    //TODO: do some fancy enable-if stuff to avoid copying DofInt to a PetscInt array if unnecessary

    void finishMatrixAssembly_impl() override
    {
      MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
    }


    // factor the matrix/update the preconditioner.
    void factor_impl() override
    {
      if (m_new_matrix_values)
      {
        std::cout << "factoring preconditioner" << std::endl;
        PCSetReusePreconditioner(m_pc, PETSC_FALSE);
        PCSetUp(m_pc);
        PCSetReusePreconditioner(m_pc, PETSC_TRUE);
        std::cout << "finished factoring preconditioner" << std::endl;
      }
    }

    // solve Ax = b for x.  The initial value of x might be used for an initial guess.
    // factor() must have been called first
    void solve_impl(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x) override
    {
      if (m_new_matrix_values)
      {
        std::cout << "about to call KSPSetUp" << std::endl;
        KSPSetUp(m_ksp);
        std::cout << "finished calling KSPSetUp" << std::endl;
        m_new_matrix_values = false;
      }

      copyVec(b, m_b);
      copyVec(x, m_x);

      KSPSolve(m_ksp, m_b, m_x);
      
      KSPConvergedReason reason;
      KSPGetConvergedReason(m_ksp, &reason);

      if (reason < 0)
      {
        const char* reason_str;
        KSPGetConvergedReasonString(m_ksp, &reason_str);
        throw std::runtime_error(std::string("KSP failed to converge for reason: ") + reason_str);
      }
      
      copyVec(m_x, x);
    }

    void matVec_impl(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      copyVec(x, m_x);
      MatMult(m_A, m_x, m_b);
      copyVec(m_b, b);
    }

  private:

    void copyVec(const ArrayType<Real, 1>& vec_in, Vec vec_out)
    {
      PetscScalar* vec_p;
      int len = vec_in.shape()[0];
      VecGetArray(vec_out, &vec_p);
      std::copy(&(vec_in[0]), &(vec_in[0]) + len, vec_p);
      VecRestoreArray(vec_out, &vec_p);
    }


    void copyVec(const Vec vec_in, ArrayType<Real, 1>& vec_out)
    {
      PetscScalar* vec_p;
      int len = vec_out.shape()[0];
      VecGetArray(vec_in, &vec_p);
      std::copy(vec_p, vec_p + len, &(vec_out[0]));
      VecRestoreArray(vec_in, &vec_p);
    }


    LargeMatrixOptsPetsc m_opts;
    Vec m_x;
    Vec m_b;
    Mat m_A;
    KSP m_ksp;
    PC m_pc;
    bool m_new_matrix_values = true;
};


} // namespace

#endif