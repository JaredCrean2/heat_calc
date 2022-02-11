#ifndef REFERENCE_ELEMENT_NODES_H
#define REFERENCE_ELEMENT_NODES_H

#include <vector>
#include "reference_element_geometry.h"

namespace reference_element {

struct ReferenceNode
{
  Point xi;
};


class ReferenceNodes
{
  public:
    using NodeVec = std::vector<ReferenceNode>;
    ReferenceNodes(const NodeVec& vert_nodes, const NodeVec& edge_nodes,
                   const NodeVec& face_nodes, const NodeVec& element_nodes) :
      m_nodes{vert_nodes, edge_nodes, face_nodes, element_nodes}
    {}

  int getNumNodes(int dim) { return m_nodes.at(dim).size(); }

  const ReferenceNode& getNode(int dim, int node) { return m_nodes.at(dim).at(node); }

  private:
    std::array<NodeVec, 4> m_nodes;
};


ReferenceNodes getLagrangeHexNodes(int degree);

}
#endif