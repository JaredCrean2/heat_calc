#ifndef REFERENCE_ELEMENT_GEOMETRY_DEF_H
#define REFERENCE_ELEMENT_GEOMETRY_DEF_H

#include <vector>
#include "mesh/reference_element_geometry.h"

namespace reference_element {

struct ReferenceElementDef
{
  using Int = int_least8_t;
  ReferenceElementDef(Int nverts, const std::vector<std::array<Int, 2>>&edge_verts,
                      const std::vector<std::vector<Int>>& face_edges,
                      const std::vector<std::vector<Int>> face_edge_perms,
                      const std::vector<Int>& el_faces,
                      const std::vector<Int>& el_face_perms,
                      const std::vector<Int>& verts_to_apf,
                      const std::vector<Int>& edges_to_apf,
                      const std::vector<Int>& faces_to_apf) :
    nverts(nverts),
    edge_verts(edge_verts),
    face_edges(face_edges),
    face_edge_perms(face_edge_perms),
    el_faces(el_faces),
    el_face_perms(el_face_perms),
    verts_to_apf(verts_to_apf),
    edges_to_apf(edges_to_apf),
    faces_to_apf(faces_to_apf)
  {}

  Int nverts;
  std::vector<std::array<Int, 2>> edge_verts;
  std::vector<std::vector<Int>> face_edges;
  std::vector<std::vector<Int>> face_edge_perms;
  std::vector<Int> el_faces;
  std::vector<Int> el_face_perms;

  // nodemap to apf refence element
  // apf only supports up to 2nd order (1 node per entity)
  // so we don't store permutation information
  std::vector<Int> verts_to_apf;
  std::vector<Int> edges_to_apf;
  std::vector<Int> faces_to_apf;
};

void checkReferenceElementDef(const ReferenceElementDef& def);

class ReferenceElementGeometry
{
  public:
    explicit ReferenceElementGeometry(const ReferenceElementDef& def);

    int getNumVerts() const { return m_verts.size(); }

    GEPtr getVert(int i) { return m_verts.at(i); }

    int getNumEdges() const { return m_edges.size(); }

    GEPtr getEdge(int i) { return m_edges.at(i); }

    int getNumFaces() const { return m_faces.size(); }

    GEPtr getFace(int i) { return m_faces.at(i); }

    GEPtr getElement() { return m_element; }

    const ReferenceElementDef& getDef() const {return m_def; }

  private:
    std::vector<GEPtr> m_verts;
    std::vector<GEPtr> m_edges;
    std::vector<GEPtr> m_faces;
    GEPtr m_element;
    ReferenceElementDef m_def;
};

}  // namespace

#endif