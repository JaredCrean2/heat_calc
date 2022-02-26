#include "mesh/reference_element_nodes.h"

namespace reference_element {

ReferenceNodes::NodeVec getTPNodes1D(const std::vector<Real>& node_xi )
{
  ReferenceNodes::NodeVec nodevec;
  for (unsigned int i=0; i < node_xi.size(); ++i)
    nodevec.push_back( ReferenceNode{ {node_xi[i], 0, 0} });

  return nodevec;
}


ReferenceNodes::NodeVec getTPNodes2D(const std::vector<Real>& node_xi )
{
  ReferenceNodes::NodeVec nodevec;
  for (unsigned int j=0; j < node_xi.size(); ++j)
    for (unsigned int i=0; i < node_xi.size(); ++i)
      nodevec.push_back( ReferenceNode{ {node_xi[i], node_xi[j], 0} });

  return nodevec;
}


ReferenceNodes::NodeVec getTPNodes3D(const std::vector<Real>& node_xi )
{
  ReferenceNodes::NodeVec nodevec;
  for (unsigned int k=0; k < node_xi.size(); ++k)
    for (unsigned int j=0; j < node_xi.size(); ++j)
      for (unsigned int i=0; i < node_xi.size(); ++i)
        nodevec.push_back( ReferenceNode{ {node_xi[i], node_xi[j], node_xi[k]} });

  return nodevec;
}


ReferenceNodes getLagrangeHexNodes(int degree)
{
  using NodeVec = ReferenceNodes::NodeVec;
  
  switch (degree)
  {
    case 1:
    {
      ReferenceNode vertnode{ {0, 0, 0} };
      NodeVec vert_nodes{vertnode}, edge_nodes, face_nodes, element_nodes;
      return ReferenceNodes(vert_nodes, edge_nodes, face_nodes, element_nodes);
    }

    case 2:
    {
      ReferenceNode vertnode{ {0, 0, 0} }, edgenode{ {0.5, 0, 0} },
                    facenode{ {0.5, 0.5, 0} }, elnode{ {0.5, 0.5, 0.5} };
      NodeVec vert_nodes{vertnode}, edge_nodes{edgenode}, face_nodes{facenode}, element_nodes{elnode};
      return ReferenceNodes(vert_nodes, edge_nodes, face_nodes, element_nodes);
    }

    case 3:
    {
      ReferenceNode vertnode{ {0, 0, 0} };
      NodeVec vert_nodes{vertnode};

      Real ot = 1.0/3.0, tt = 2.0/3.0;
      std::vector<Real> node_xi{ot, tt};
      auto edge_nodes    = getTPNodes1D(node_xi);
      auto face_nodes    = getTPNodes2D(node_xi);
      auto element_nodes = getTPNodes3D(node_xi);
                    
      return ReferenceNodes(vert_nodes, edge_nodes, face_nodes, element_nodes);
    }

    default:
      throw std::invalid_argument(std::string("unsupported degree in getLagrangeHex(): ") + std::to_string(degree));
  };
}

}