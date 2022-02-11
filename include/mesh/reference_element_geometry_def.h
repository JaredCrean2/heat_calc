#ifndef REFERENCE_ELEMENT_GEOMETRY_DEF_H
#define REFERENCE_ELEMENT_GEOMETRY_DEF_H

#include <vector>
#include "mesh/reference_element_geometry.h"

namespace reference_element {

struct ReferenceElementDef
{
  using Int = int_least8_t;
  Int nverts;
  std::vector<std::array<Int, 2>> edge_verts;
  std::vector<std::vector<Int>> face_edges;
  std::vector<std::vector<Int>> face_edge_perms;
  std::vector<Int> el_faces;
  std::vector<Int> el_face_perms;
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

  private:
    std::vector<GEPtr> m_verts;
    std::vector<GEPtr> m_edges;
    std::vector<GEPtr> m_faces;
    GEPtr m_element;
};

}  // namespace

#endif