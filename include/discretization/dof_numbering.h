#ifndef DISC_DOF_NUMBERING_H
#define DISC_DOF_NUMBERING_H

#include "mesh/mesh.h"
#include "discretization/volume_discretization.h"

struct ElementNode
{
  ElementNode(Index element, LocalIndex node) :
    element(element),
    node(node)
{}

  Index element;    // local element number
  LocalIndex node;  // local node number
};

class DofNumbering
{
  public:
    explicit DofNumbering(std::shared_ptr<Mesh::MeshCG> mesh) :
    m_num_dofs(mesh->getNumDofs())
    {
      for (int i=0; i < mesh->getNumVolumeGroups(); ++i)
        setDofs(mesh, i);
    }


    const ArrayType<Index, 2>& getDofs(const int idx) const
    {
      return m_dof_nums.at(idx);
    }

    const ArrayType<Index, 2>& getDofs(const VolDiscPtr vol_disc) const
    {
      return getDofs(vol_disc->getIdx());
    }

    const std::vector<ElementNode>& getDirichletNodes(const int idx) const
    {
      return m_dirichlet_node_nums.at(idx);
    }

    const std::vector<ElementNode>& getDirichletNodes(const VolDiscPtr vol_disc) const
    {
      return getDirichletNodes(vol_disc->getIdx());
    }


    bool isDofActive(Index dof) const
    {
      return dof < m_num_dofs;
    }

  private:
    
    void setDofs(std::shared_ptr<Mesh::MeshCG> mesh, int vol_idx);

    std::vector<ArrayType<Index, 2>> m_dof_nums;
    std::vector<std::vector<ElementNode>> m_dirichlet_node_nums;
    int m_num_dofs;

};

using DofNumberingPtr = std::shared_ptr<DofNumbering>;


#endif