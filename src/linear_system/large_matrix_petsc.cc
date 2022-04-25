#include "linear_system/large_matrix_petsc.h"

namespace linear_system {

LargeMatrixPetsc::LargeMatrixPetsc(DofInt mlocal, DofInt nlocal, LargeMatrixOptsPetsc opts, std::shared_ptr<SparsityPattern> sparsity_pattern) :
  LargeMatrix(mlocal, nlocal),
  m_opts(opts)
{  
  bool is_symmetric_mattype = false;
  if (opts.petsc_opts.count("mat_type") == 0)
  {
    std::cout << "setting mat type" << std::endl;
    opts.petsc_opts["mat_type"] = opts.is_value_symmetric ? "sbaij" : "aij";
  } else
  {
    if (opts.petsc_opts["mat_type"] == "sbaij")
      is_symmetric_mattype = true;
  }  
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
  //TODO: what happens if MAT_SYMMETRIC is true but a non-symmetric matrix format is used?
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

  const PetscInt *diagonal_counts = nullptr, *offproc_counts=nullptr;
  const PetscInt *diagonal_counts_sym = nullptr, *offproc_counts_sym=nullptr;
  if (is_symmetric_mattype)
  {
    diagonal_counts_sym = sparsity_pattern->getDiagonalCountsSym().data();
    offproc_counts_sym  = sparsity_pattern->getOffProcCountsSym().data();
  } else
  {
    diagonal_counts = sparsity_pattern->getDiagonalCounts().data();
    offproc_counts  = sparsity_pattern->getOffProcCounts().data();
  }
    
  MatXAIJSetPreallocation(m_A, 1, diagonal_counts, offproc_counts,
                                  diagonal_counts_sym, offproc_counts_sym);
  

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


LargeMatrixPetsc::~LargeMatrixPetsc()
{
  VecDestroy(&m_x);
  VecDestroy(&m_b);
  MatDestroy(&m_A);
  KSPDestroy(&m_ksp);
  // KSP manages the PC, no need to destroy it explicitly
}


void LargeMatrixPetsc::factor_impl()
{
  if (m_new_matrix_values)
  {
    PCSetReusePreconditioner(m_pc, PETSC_FALSE);
    PCSetUp(m_pc);
    PCSetReusePreconditioner(m_pc, PETSC_TRUE);
  }
}


void LargeMatrixPetsc::solve_impl(const ArrayType<Real, 1>& b, ArrayType<Real, 1>& x)
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

void LargeMatrixPetsc::matVec_impl(const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  copyVec(x, m_x);
  MatMult(m_A, m_x, m_b);
  copyVec(m_b, b);
}


void LargeMatrixPetsc::copyVec(const ArrayType<Real, 1>& vec_in, Vec vec_out)
{
  PetscScalar* vec_p;
  int len = vec_in.shape()[0];
  VecGetArray(vec_out, &vec_p);
  std::copy(&(vec_in[0]), &(vec_in[0]) + len, vec_p);
  VecRestoreArray(vec_out, &vec_p);
}


void LargeMatrixPetsc::copyVec(const Vec vec_in, ArrayType<Real, 1>& vec_out)
{
  PetscScalar* vec_p;
  int len = vec_out.shape()[0];
  VecGetArray(vec_in, &vec_p);
  std::copy(vec_p, vec_p + len, &(vec_out[0]));
  VecRestoreArray(vec_in, &vec_p);
}


void setPetscOptions(const LargeMatrixOptsPetsc& opts)
{
  for (auto& p : opts.petsc_opts)
  {
    std::string full_key = std::string("-") + opts.opts_prefix + p.first;
    PetscOptionsSetValue(NULL, full_key.c_str(), p.second.c_str());
  }
}

}