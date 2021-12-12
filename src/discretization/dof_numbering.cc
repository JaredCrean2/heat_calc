#include <cassert>
#include "discretization/dof_numbering.h"

void DofNumbering::setDofs(std::shared_ptr<Mesh::MeshCG> mesh, int vol_idx)
{
  assert(m_dof_nums.size() < vol_idx);

  auto& vol_group = mesh->getElements(vol_idx);
  m_dof_nums.emplace_back(boost::extents[vol_group.getNumElems()][vol_group.getNumSolPtsPerElement()]);
  m_dirichlet_node_nums.emplace_back();

  auto& dof_nums = m_dof_nums.back();
  auto& dirichlet_dof_nums = m_dirichlet_node_nums.back();
  std::vector<Index> nodenums_tmp;
  std::set<Index> dirchlet_dofs_seen;

  for (int i=0; i < vol_group.getNumElems(); ++i)
  {
    mesh->getElementDofs(vol_group, i, nodenums_tmp);

    assert(nodenums_tmp.size() == vol_group.getNumSolPtsPerElement());

    for (int j=0; j < vol_group.getNumSolPtsPerElement(); ++j)
    {
      auto dof = nodenums_tmp[j];
      dof_nums[i][j] = dof;

      if ( !(mesh->isDofActive(dof)) && dirchlet_dofs_seen.count(dof) == 0 )  // if dirichlet
      {
        dirichlet_dof_nums.emplace_back(i, j);
        dirchlet_dofs_seen.insert(dof);
      }
    }
  }
}