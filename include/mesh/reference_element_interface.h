#ifndef REFERENCE_ELEMENT_DEF_H
#define REFERENCE_ELEMENT_DEF_H

#include "mesh/reference_element_geometry.h"
#include "mesh/reference_element_geometry_interface.h"
//#include "mesh/reference_element_geometry_hex.h"
#include "mesh/reference_element_nodes.h"
#include "utils/math.h"

namespace reference_element {

struct NodeLocation
{
  int dim;
  int entity_index;
  int node_index;
};

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

    // get total number of nodes on all entities <= dim
    int getNumNodesInclusive(int dim)
    {
      int total = 0;
      for (int i=0; i <= dim; ++i)
        total += getNumEntities(i) * getNumNodes(i);

      return total;
    }

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

    // given a node in the range [0, getNumNodesTotal), returns the
    // {dimension, entity index within dimension, node on entity}
    // entity is in the RE ordering
    NodeLocation getNodeLocation(int node)
    {
      int offset = 0, dim;
      for (dim=0; dim <= 3; ++dim)
      {
        int offset_tmp = getNumEntities(dim) * getNumNodes(dim);
        if (node < offset + offset_tmp)
          break;

        offset += offset_tmp;
      }

      int node_on_dim = node - offset;
      int entity = node_on_dim / getNumNodes(dim);
      int node_on_entity = node_on_dim % getNumNodes(dim);

      return {dim, entity, node_on_entity};
    }

    // given the dimension and apf index, get the ReferenceElement index
    int getREEntityIndex(int dim, int idx);
    
    //TODO: this duplicated getApfEntity above?
    // given the dimension and ReferenceElement index, get the apf entity index
    int getApfEntityIndex(int dim, int idx)
    {
      switch (dim)
      {
        case 0: return getDef().verts_to_apf.at(idx);
        case 1: return getDef().edges_to_apf.at(idx);
        case 2: return getDef().faces_to_apf.at(idx);
        case 3: return idx;
        default:
          throw std::runtime_error(std::string("Unrecognized dimension: ") + std::to_string(dim));
      }
    }

    // returns nfaces x nnodes per face array giving the indices of the nodes on each face
    const ArrayType<LocalIndex, 2>& getFaceNodes() const { return m_face_nodes; }

    int getNumFaceNodes() const { return getFaceNodes().shape()[1]; }

    int getNumNodesTP() const { return m_tp_nodemap.shape()[0]; }

    int getDegree() const { return getNumNodesTP() - 1; }

    // getNumNodesTotal x 3
    const ArrayType<Real, 2>& getNodeXi() const { return m_node_xi; }

    // nfaces x 3
    const ArrayType<Real, 2>& getNormals() const { return m_normals; }

    // 3D array, each dimension getNumNodesTP
    const ArrayType<LocalIndex, 3>& getTPNodemap() const {return m_tp_nodemap; }

    // length getNumNodesT() 
    const std::vector<Real>& getTensorProductXi() const { return m_tp_xi; }

    std::pair<Real, Real> getXiRange() const { return std::pair<Real, Real>(m_tp_xi.front(), m_tp_xi.back()); }


  private:

    ReferenceElementGeometry m_geom;
    ReferenceNodes m_nodes;
    ArrayType<LocalIndex, 2> m_face_nodes;
    ArrayType<LocalIndex, 3> m_tp_nodemap;
    ArrayType<Real, 2> m_normals;
    std::vector<Real> m_tp_xi;
    ArrayType<Real, 2> m_node_xi;  // num_total_nodes x 3
};

namespace Impl {

void getFaceNodes(ReferenceElement& ref_el, ArrayType<LocalIndex, 2>& face_nodes);

void getNodeXi(ReferenceElement& ref_el, ArrayType<Real, 2>& node_xi);

std::vector<Real> get1DXi(ReferenceElement& ref_el);

int searchNode(ReferenceElement& ref_el, const Point& pt);

void computeTPNodemap(ReferenceElement& ref_el, ArrayType<LocalIndex, 3>& nodemap);

void computeNormals(ReferenceElement& ref_el, ArrayType<Real, 2>& normals);
}

using REPtr = std::shared_ptr<ReferenceElement>;

REPtr getLagrangeHexReferenceElement(int degree);


}

#endif
