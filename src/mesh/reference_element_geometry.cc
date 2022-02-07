#include "mesh/reference_element_geometry.h"
#include "mesh/reference_element_geometry_types.h"
#include <iostream>

namespace reference_element {

void GeometricEntity::setDown(GEPtr entity, int perm)
{ 
  assertAlways(entity->getDimension() == getDimension() - 1, "Downward adjacency must be one dimension lower");
  m_down.push_back(entity);
  m_down_perm.push_back(perm);
}


Point GeometricEntity::computeFromDownCoordinates(int i, const Point pt)
{
  assertAlways(i >= 0 && static_cast<unsigned int>(i) <= m_down.size(), "invalid downward entity");
  validatePoint(pt);
  return computefromDownCoordinatesImpl(i, pt);
}


Point GeometricEntity::computeFromDownCoordinates(GEPtr entity, const Point pt)
{
  for (size_t i=0; i < m_down.size(); ++i)
    if (entity == m_down[i])
      return computeFromDownCoordinates(i, pt);

  throw std::invalid_argument("entity must be a downward adjacent entity");
}


void GeometricEntity::validatePoint(const Point pt)
{
  const Real tol = 1e-13;
  for (int i=0; i < 3; ++i)
  {
    if (i >= getDimension())
      assertAlways(std::abs(pt[i]) < tol, "extra dimensions must have zero coordinate");
    else
      assertAlways(pt[i] >= -tol && pt[i] <= 1 + tol, "coordinate is out of range");
  }
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


ReferenceElementDef getStandardReferenceElementDef()
{
  return ReferenceElementDef();
}


// returns true if to_entity is a downward adjacency of from_entity 
// (traversing several dimensions down if necessary)
bool hasDownwardAdjacency(GEPtr from_entity, GEPtr to_entity)
{
  assertAlways(from_entity->getDimension() > to_entity->getDimension(), 
          "to_entity must be lower dimension than from_entity");
  
  // recursive base case
  if (from_entity->getDimension() - to_entity->getDimension() == 1)
  {
    for (int i=0; i < from_entity->countDown(); ++i)
    {
      if (from_entity->getDown(i) == to_entity)
        return true;
    }

    return false;
  } else
  {
    for (int i=0; i < from_entity->countDown(); ++i)
      if (hasDownwardAdjacency(from_entity->getDown(i), to_entity))
        return true;

    return false;
  }
}


// T must be an array types
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



//TODO: this is very slow.  See if there is a way to combine hasDownwardAdjacency
//      into this function, rather than evaluating it so many times
// gets sequence of entities that connect from_entity to to_entity such that
// entity i+1 is on the closure of entity i
// The returned sequence is in order of decreasing dimension, and from_entity
// does is the first item in the list.  to_entity does not appear in the sequence,
// the final entity is of one dimension higher.
std::vector<GEPtr> getEntityChain(GEPtr from_entity, GEPtr to_entity)
{
  assertAlways(hasDownwardAdjacency(from_entity, to_entity), 
           "cannot find chain of entities connecting given entities");

  std::vector<GEPtr> entities;
  
  GEPtr current_entity = from_entity;
  entities.push_back(current_entity);
  while (current_entity->getDimension() - to_entity->getDimension() > 1)
  {
    bool found = false;
    for (int i=0; i < current_entity->countDown(); ++i)
    {
      if (hasDownwardAdjacency(current_entity->getDown(i), to_entity))
      {
        current_entity = current_entity->getDown(i);
        entities.push_back(current_entity);
        found = true;
        break;
      }
    }

    // strictly speaking, this check is unnecessary because the check at the beginning of
    // this function guarantees hasDownwardAdjacency() will return true for some input value
    if (!found)
      throw std::runtime_error("cannot find chain of entities connecting given entities");
  }

  return entities;
}


Point reclassifyPoint(GEPtr from_entity, const Point& pt_in, GEPtr to_entity)
{
  assertAlways(from_entity->getDimension() < to_entity->getDimension(), 
          "can only reclassify point from lower dimension entity to higher dimension entity");

  auto entities = getEntityChain(to_entity, from_entity);

  Point pt = pt_in;
  GEPtr prev_entity = from_entity;
  for (auto it = entities.rbegin(); it != entities.rend(); ++it)
  {
    pt = (*it)->computeFromDownCoordinates(prev_entity, pt);
    prev_entity = *it;
  }

  return pt;
}


} // namespace