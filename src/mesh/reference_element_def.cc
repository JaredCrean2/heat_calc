#include "mesh/reference_element_def.h"
#include "mesh/reference_element_geometry_types.h"

namespace reference_element {

// T must be an array type
// This checks that no elements are permutive equivalent
template <typename T>
bool checkUnique(std::vector<T> vals)
{
  for (auto& v : vals)
    std::sort(v.begin(), v.end());

  return std::unique(vals.begin(), vals.end()) == vals.end();
}

void checkReferenceElementDef(const ReferenceElementDef& def)
{
  // check uniquenesss
  assertAlways(checkUnique(def.edge_verts), "Edge definitions are not unique");
  assertAlways(checkUnique(def.face_edges), "Face definitions are not unique");

  auto el_faces_copy = def.el_faces;
  assertAlways(std::unique(el_faces_copy.begin(), el_faces_copy.end()) == el_faces_copy.end(),
               "Element definition not unique"); 

  // check sizes
  assertAlways(def.face_edge_perms.size() == def.face_edges.size(), "face edges permutation size is incorrect");
  assertAlways(def.el_face_perms.size() == def.el_faces.size(), "face edges permutation size is incorrect");

  // check every entity is used
  std::vector<bool> vert_used(def.nverts), edge_used(def.edge_verts.size()), face_used(def.face_edges.size());
  for (auto& verts : def.edge_verts)
  {
    vert_used.at(verts[0]) = true;
    vert_used.at(verts[1]) = true;
  }

  for (auto& edges : def.face_edges)
    for (size_t j=0; j < edges.size(); ++j)
      edge_used.at(edges[j]) = true;

  for (auto& face : def.el_faces)
    face_used.at(face) = true;

  auto pred = [](const bool& b) { return b; };
  assertAlways(std::all_of(vert_used.begin(), vert_used.end(), pred),
               "Not all vertices are used in edge definitions");
  assertAlways(std::all_of(edge_used.begin(), edge_used.end(), pred),
               "Not all edges are used in face definitions");
  assertAlways(std::all_of(face_used.begin(), face_used.end(), pred),
               "Not all faces are used in element definition");
}

ReferenceElement::ReferenceElement(const ReferenceElementDef& def)
{
  checkReferenceElementDef(def);

  for (int i=0; i < def.nverts; ++i)
    m_verts.push_back(std::make_shared<GeometricVertex>(0));

  for (auto& verts_i : def.edge_verts)
  {
    auto entity = std::make_shared<GeometricEdge>(1);
    entity->setDown(m_verts[verts_i[0]]);
    entity->setDown(m_verts[verts_i[1]]);
    m_edges.push_back(entity);
  }

  for (size_t i=0; i < def.face_edges.size(); ++i)
  {
    auto entity = std::make_shared<GeometricQuad>(2);
    auto& edge_idxs = def.face_edges[i];
    auto& edge_perms = def.face_edge_perms[i];

    for (size_t j=0; j < edge_idxs.size(); ++j)
      entity->setDown(m_edges[edge_idxs[j]], edge_perms[j]);

    m_faces.push_back(entity);
  }

  m_element = std::make_shared<GeometricHex>(3);
  for (size_t j=0; j < def.el_faces.size(); ++j)
    m_element->setDown(m_faces[def.el_faces[j]], def.el_face_perms[j]);
}

}