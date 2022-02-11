#ifndef REFERENCE_ELEMENT_DEF_H
#define REFERENCE_ELEMENT_DEF_H

#include "mesh/reference_element_geometry_def.h"
#include "mesh/reference_element_nodes.h"

namespace reference_element {

class ReferenceElement
{
  public:
    ReferenceElement(ReferenceElementGeometry geom, ReferenceNodes nodes) :
      m_geom(geom),
      m_nodes(nodes)
    {}

    int getNumVerts() const { return m_geom.getNumVerts(); }

    GEPtr getVert(int i) { return m_geom.getVert(i); }

    int getNumEdges() const { return m_geom.getNumEdges(); }

    GEPtr getEdge(int i) { return m_geom.getEdge(i); }

    int getNumFaces() const { return m_geom.getNumFaces(); }

    GEPtr getFace(int i) { return m_geom.getFace(i); }

    GEPtr getElement() { return m_geom.getElement(); }

    int getNumEntities(int dim)
    {
      switch (dim)
      {
        case 0: return getNumVerts();
        case 1: return getNumEdges();
        case 2: return getNumFaces();
        case 3: return 1;
        default:
          throw std::invalid_argument("invalid entity dimension");
      }
    }

    // gets number of nodes on entities of given dimension
    int getNumNodes(int dim) { return m_nodes.getNumNodes(dim); }

    // total number of nodes on all entities
    int getNumNodesTotal() { return getNumNodesInclusive(3); }

    // get the index of the node in the range [0, getNumNodesTotal)
    // node is in the range [0, total number of nodes on all entities of given dim)
    int getNodeIndex(int dim, int node)
    {
      assert(node >= 0 && node < getNumNodes(dim) * getNumEntities(dim));
      return getNumNodesInclusive(dim - 1) + node;
    }

    int getNodeIndex(int dim, int entity, int node)
    {
      assert(entity >= 0 && entity < getNumEntities(dim));
      assert(node >= 0   && node < getNumNodes(dim));

      return getNodeIndex(dim, entity * getNumNodes(dim) + node);
    }

  private:
    // get total number of nodes on all entities <= dim
    int getNumNodesInclusive(int dim)
    {
      int total = 0;
      for (int i=0; i <= dim; ++i)
        total += getNumEntities(i) * getNumNodes(i);

      return total;
    }

    ReferenceElementGeometry m_geom;
    ReferenceNodes m_nodes;
};

}

#endif