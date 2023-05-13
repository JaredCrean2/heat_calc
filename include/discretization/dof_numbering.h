#ifndef DISC_DOF_NUMBERING_H
#define DISC_DOF_NUMBERING_H

#include "discretization/volume_discretization.h"
#include "mesh/mesh.h"

namespace Mesh {
  class MeshCG;
}

struct ElementNode
{
  ElementNode(Index element, LocalIndex node) :
    element(element),
    node(node)
{}

  Index element;    // local element number
  LocalIndex node;  // local node number
};

template <typename T>
struct Assign2
{
  T operator()(const T first, const T second) { return second;}
};

class DofNumbering
{
  public:
    explicit DofNumbering(std::shared_ptr<Mesh::MeshCG> mesh);

    // returns array numElems in volume group x numSolPtsPerElement giving
    // dof numbers for given block
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
      return dof < m_num_local_dofs;
    }

    int getNumLocalDofs() const
    {
      return m_num_local_dofs;
    }

    int getNumOwnedDofs() const
    {
      return m_num_owned_dofs;
    }

    int getNumDirichletNodeSections() const { return m_dest_node_ranges.size() - 1;}

    Mesh::NodeTriplet getSrcDirichletNode(int section) const { return m_src_dirichlet_nodes[section]; }

    void getDestDirichletNodes(int section, std::vector<Mesh::NodeTriplet>& dest_nodes)
    {
      int start_idx = m_dest_node_ranges[section];
      int end_idx = m_dest_node_ranges[section+1];
      dest_nodes.resize(end_idx - start_idx);
      for (int i=start_idx; i < end_idx; ++i)
        dest_nodes[i - start_idx] = m_dest_dirichlet_nodes[i];
    }

  private:
    
    void setDofs(std::shared_ptr<Mesh::MeshCG> mesh, int vol_idx);

    std::vector<ArrayType<Index, 2>> m_dof_nums;
    std::vector<std::vector<ElementNode>> m_dirichlet_node_nums;  //TODO: unused?
    std::vector<Mesh::NodeTriplet> m_src_dirichlet_nodes;
    std::vector<Index> m_dest_node_ranges;
    std::vector<Mesh::NodeTriplet> m_dest_dirichlet_nodes;
    int m_num_local_dofs;
    int m_num_owned_dofs;

};

using DofNumberingPtr = std::shared_ptr<DofNumbering>;


#endif