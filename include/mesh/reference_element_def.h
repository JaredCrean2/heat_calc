#ifndef REFERENCE_ELEMENT_DEF_H
#define REFERENCE_ELEMENT_DEF_H

#include "mesh/reference_element_geometry.h"
#include "mesh/reference_element_geometry_def.h"
#include "mesh/reference_element_geometry_types.h"
#include "mesh/reference_element_nodes.h"
#include "utils/math.h"

namespace reference_element {

class ReferenceElement
{
  public:
    ReferenceElement(const ReferenceElementGeometry& geom, const ReferenceNodes& nodes);

    const ReferenceElementDef& getDef() const {return m_geom.getDef(); }

    int getNumVerts() const { return m_geom.getNumVerts(); }

    GEPtr getVert(int i) { return m_geom.getVert(i); }

    int getNumEdges() const { return m_geom.getNumEdges(); }

    GEPtr getEdge(int i) { return m_geom.getEdge(i); }

    int getNumFaces() const { return m_geom.getNumFaces(); }

    GEPtr getFace(int i) { return m_geom.getFace(i); }

    GEPtr getElement() { return m_geom.getElement(); }

    GEPtr getEntity(int dim, int idx)
    {
      switch (dim)
      {
        case 0: return getVert(idx);
        case 1: return getEdge(idx);
        case 2: return getFace(idx);
        case 3: return getElement();
        default:
          throw std::invalid_argument("invalid entity dimension");
      }
    }


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

    int getApfEntity(int dim, int entity)
    {
      switch (dim)
      {
        case 0: return m_geom.getDef().verts_to_apf.at(entity);
        case 1: return m_geom.getDef().edges_to_apf.at(entity);
        case 2: return m_geom.getDef().faces_to_apf.at(entity);
        case 3: return 0;
        default:
          throw std::invalid_argument("invalid entity dimension");
      }  
    }

    // gets number of nodes on entities of given dimension
    int getNumNodes(int dim) { return m_nodes.getNumNodes(dim); }

    ReferenceNode getNode(int dim, int node) { return m_nodes.getNode(dim, node); }

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

    const ArrayType<LocalIndex, 3>& getTPNodemap() const {return m_tp_nodemap; }

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
    ArrayType<LocalIndex, 3> m_tp_nodemap;
    ArrayType<Real, 2> m_normals;
};


std::vector<Real> get1DXi(ReferenceElement& ref_el);

int searchNode(ReferenceElement& ref_el, const Point& pt);

ArrayType<LocalIndex, 3> computeTPNodemap(ReferenceElement& ref_el);

ArrayType<Real, 2> computeNormals(ReferenceElement& ref_el);

std::shared_ptr<ReferenceElement> getLagrangeHexReferenceElement(int degree);

}

#endif