#ifndef SPARSITY_PATTERN_MESH_H
#define SPARSITY_PATTERN_MESH_H

#include "linear_system/sparsity_pattern.h"
#include "mesh/mesh.h"

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

    void computePattern(bool symmetric)
    {
      auto& local_dofs  = symmetric ? m_local_dofs_sym : m_local_dofs;
      auto& remote_dofs = symmetric ? m_remote_dofs_sym : m_remote_dofs_sym;

      local_dofs.resize(m_mesh->getNumDofs());
      remote_dofs.resize(m_mesh->getNumDofs());

      std::vector<DofInt> el_dofs, connected_dofs;
      for (int i=0; i < m_mesh->getNumVolumeGroups(); ++i)
      {
        //std::cout << "\nvolume group " << i << std::endl;
        auto vol_group = m_mesh->getElements(i);
        for (int el=0; el < vol_group.getNumElems(); ++el)
        {
          m_mesh->getElementDofs(vol_group, el, el_dofs);
          for (int j=0; j < vol_group.getNumSolPtsPerElement(); ++j)
            if (m_mesh->isDofActive(el_dofs[j]))
            {
              //TODO: skip dofs already processed
              //std::cout << "element " << el << ", node " << j << std::endl;
              DofInt dof_j = el_dofs[j];
              m_mesh->getDofConnectivity(vol_group, el, j, connected_dofs);

              int count = 0;
              if (symmetric)
              {
                // count dofs on diagonal + upper triangle
                for (auto& dof : connected_dofs)
                  count += dof >= dof_j && m_mesh->isDofActive(dof) ? 1 : 0;
              } else
              {
                for (auto& dof : connected_dofs)
                  count += m_mesh->isDofActive(dof) ? 1 : 0;
              }

              local_dofs[dof_j]  = count;
              remote_dofs[dof_j] = 0;
            }
        }
      }

      if (symmetric)
        m_computed_symmetric = true;
      else
        m_computed_nonsymmetric = true;
    }

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