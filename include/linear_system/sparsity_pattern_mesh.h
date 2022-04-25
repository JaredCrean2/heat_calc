#ifndef SPARSITY_PATTERN_MESH_H
#define SPARSITY_PATTERN_MESH_H

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
    SparsityPatternMesh(std::shared_ptr<Mesh::MeshCG> mesh) :
      m_mesh(mesh)
    {}

    const std::vector<PetscInt>& getDiagonalCounts() override
    { 
      if (!m_computed_nonsymmetric)
        computePattern(false);

      return m_local_dofs;
    }

    const std::vector<PetscInt>& getDiagonalCountsSym() override
    { 
      if (!m_computed_symmetric)
        computePattern(true);

      return m_local_dofs_sym;
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

  private:
    void computePattern(bool symmetric);

    bool m_computed_nonsymmetric = false;
    bool m_computed_symmetric    = false;
    std::shared_ptr<Mesh::MeshCG> m_mesh;
    std::vector<PetscInt> m_local_dofs;
    std::vector<PetscInt> m_remote_dofs;
    std::vector<PetscInt> m_local_dofs_sym;
    std::vector<PetscInt> m_remote_dofs_sym;
};

} // namespace

#endif