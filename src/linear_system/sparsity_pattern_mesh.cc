#include "linear_system/sparsity_pattern_mesh.h"

namespace linear_system {

void SparsityPatternMesh::computePattern(bool symmetric)
{
  auto& local_dofs  = symmetric ? m_local_dofs_sym : m_local_dofs;
  auto& remote_dofs = symmetric ? m_remote_dofs_sym : m_remote_dofs_sym;

  local_dofs.resize(m_mesh->getNumDofs());
  remote_dofs.resize(m_mesh->getNumDofs());
  std::fill(local_dofs.begin(), local_dofs.end(), -1);

  std::vector<DofInt> el_dofs, connected_dofs;
  for (int i=0; i < m_mesh->getNumVolumeGroups(); ++i)
  {
    //std::cout << "\nvolume group " << i << std::endl;
    auto vol_group = m_mesh->getElements(i);
    for (int el=0; el < vol_group.getNumElems(); ++el)
    {
      m_mesh->getElementDofs(vol_group, el, el_dofs);
      for (int j=0; j < vol_group.getNumSolPtsPerElement(); ++j)
        if (m_mesh->isDofActive(el_dofs[j]) && local_dofs[el_dofs[j]] == -1)
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

}  // namespace