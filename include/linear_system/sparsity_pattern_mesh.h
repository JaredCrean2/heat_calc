#ifndef SPARSITY_PATTERN_MESH_H
#define SPARSITY_PATTERN_MESH_H

#include "ProjectDefs.h"
#include "linear_system/sparsity_pattern.h"
#include <memory>

namespace Mesh
{
  class MeshCG;
}

namespace linear_system {

class SparsityPatternMesh : public SparsityPattern
{
  public:
    explicit SparsityPatternMesh(std::shared_ptr<Mesh::MeshCG> mesh);

    const std::vector<PetscInt>& getDiagonalCounts() override
    { 
      if (!m_computed_nonsymmetric)
        computePattern(false);

      return m_onproc_dofs;
    }

    const std::vector<PetscInt>& getDiagonalCountsSym() override
    { 
      if (!m_computed_symmetric)
        computePattern(true);

      return m_onproc_dofs_sym;
    }

    const std::vector<PetscInt>& getOffProcCounts() override 
    {
      if (!m_computed_nonsymmetric)
        computePattern(false);

      return m_remote_dofs; 
    }

    const std::vector<PetscInt>& getOffProcCountsSym() override 
    {
      if (!m_computed_symmetric)
        computePattern(true);

      return m_remote_dofs_sym; 
    }

    const std::vector<PetscInt>& getGhostGlobalIndices() override
    {
      return m_ghost_global_dofs;
    }

    const std::vector<PetscInt>& getGhostLocalIndices() override
    {
      return m_ghost_onproc_dofs;
    }

    const std::vector<PetscInt>& getOwnedToLocalInfo() override
    {
      return m_owned_dof_to_local;
    }


  private:
    void computePattern(bool symmetric);

    void getDofStatus(std::vector<DofInt>& local_to_owned_dof);

    bool m_computed_nonsymmetric = false;
    bool m_computed_symmetric    = false;
    std::shared_ptr<Mesh::MeshCG> m_mesh;
    std::vector<PetscInt> m_onproc_dofs;
    std::vector<PetscInt> m_remote_dofs;
    std::vector<PetscInt> m_onproc_dofs_sym;
    std::vector<PetscInt> m_remote_dofs_sym;
    std::vector<PetscInt> m_ghost_global_dofs;
    std::vector<PetscInt> m_ghost_onproc_dofs;
    std::vector<PetscInt> m_owned_dof_to_local;
};

} // namespace

#endif