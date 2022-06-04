#include "linear_system/sparsity_pattern_mesh.h"
#include "mesh/mesh.h"


namespace linear_system {

SparsityPatternMesh::SparsityPatternMesh(std::shared_ptr<Mesh::MeshCG> mesh) :
  m_mesh(mesh)
{
  mesh->getGhostDofInfo(m_ghost_global_dofs, m_ghost_onproc_dofs); 
  mesh->getOwnedLocalDofInfo(m_owned_dof_to_local);  //TODO: use this in getDofStatus?
}

void SparsityPatternMesh::computePattern(bool symmetric)
{
  auto& onproc_dofs  = symmetric ? m_onproc_dofs_sym : m_onproc_dofs;
  auto& remote_dofs = symmetric ? m_remote_dofs_sym : m_remote_dofs;

  onproc_dofs.resize(m_mesh->getNumDofs());
  remote_dofs.resize(m_mesh->getNumDofs());
  std::fill(onproc_dofs.begin(), onproc_dofs.end(), -1);

  std::vector<DofInt> local_to_owned_dof;
  getDofStatus(local_to_owned_dof);


  std::vector<DofInt> el_dofs, connected_dofs;
  for (int i=0; i < m_mesh->getNumVolumeGroups(); ++i)
  {
    auto vol_group = m_mesh->getElements(i);
    for (int el=0; el < vol_group.getNumElems(); ++el)
    {
      m_mesh->getElementDofs(vol_group, el, el_dofs);      
      for (int j=0; j < vol_group.getNumSolPtsPerElement(); ++j)
      {
        int local_dof_num = el_dofs[j];
        if (!m_mesh->isDofActive(local_dof_num))
          continue;

        int owned_dof_num = local_to_owned_dof[local_dof_num];
        if (m_mesh->isDofActive(el_dofs[j]) && owned_dof_num != -1 && onproc_dofs[owned_dof_num] == -1)
        {
          m_mesh->getDofConnectivity(vol_group, el, j, connected_dofs);

          int onproc_count = 0, offproc_count=0;
          if (symmetric)
          {
            // count dofs on diagonal + upper triangle
            for (auto& dof : connected_dofs)
              if (m_mesh->isDofActive(dof) && dof >= local_dof_num)
                local_to_owned_dof[dof] == -1 ? offproc_count++ : onproc_count++;
          } else
          {
            for (auto& dof : connected_dofs)
              if (m_mesh->isDofActive(dof))
                local_to_owned_dof[dof] == -1 ? offproc_count++ : onproc_count++;
          }

          onproc_dofs[owned_dof_num] = onproc_count;
          remote_dofs[owned_dof_num] = offproc_count;
        }
      }
    }
  }

  if (symmetric)
    m_computed_symmetric = true;
  else
    m_computed_nonsymmetric = true;
}

// gives a vector v such that v[local_dof] = owned_dof, or -1 if dof is not owned
void SparsityPatternMesh::getDofStatus(std::vector<DofInt>& local_to_owned_dof)
{
  local_to_owned_dof.resize(m_mesh->getNumDofs(), -1);
  std::vector<DofInt> owned_to_local_dofs;
  m_mesh->getOwnedLocalDofInfo(owned_to_local_dofs);  // TODO: maybe the mesh should cache this and return a reference
  for (size_t i=0; i < owned_to_local_dofs.size(); ++i)
    local_to_owned_dof[owned_to_local_dofs[i]] = i;
}

}  // namespace