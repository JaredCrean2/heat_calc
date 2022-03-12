#include "mesh/reference_element_geometry.h"
#include "mesh/reference_element_geometry_hex.h"
#include <iostream>
#include <set>

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


std::vector<GEPtr> getDownward(GEPtr entity, int dim)
{
  assertAlways(dim < entity->getDimension(), "invalid dimension");

  std::vector<GEPtr> down(entity->countDown());
  for (int i=0; i < entity->countDown(); ++i)
    down[i] = entity->getDown(i);
    
  if (dim == entity->getDimension() - 1)
  {
    return down;
  } else
  {
    std::vector<GEPtr> down_total;
    std::set<GEPtr> seen_entities;
    for (int i=0; i < entity->countDown(); ++i)
    {
      auto down_i = getDownward(down[i], dim);
      for (unsigned int j=0; j < down_i.size(); ++j)
        if (seen_entities.count(down_i[j]) == 0)
        {
          down_total.push_back(down_i[j]);
          seen_entities.insert(down_i[j]);
        }
    }

    return down_total;
  }
}


bool getEntityChain(GEPtr from_entity, GEPtr to_entity, std::vector<GEPtr>& entities)
{
  assertAlways(from_entity->getDimension() > to_entity->getDimension(), 
          "to_entity must be lower dimension than from_entity");

  // recursive base case
  if (from_entity->getDimension() - to_entity->getDimension() == 1)
  {
    for (int i=0; i < from_entity->countDown(); ++i)
    {
      if (from_entity->getDown(i) == to_entity)
      {
        entities.push_back(from_entity);
        return true;
      }
    }

    return false;
  } else
  {
    for (int i=0; i < from_entity->countDown(); ++i)
    {
      bool found = getEntityChain(from_entity->getDown(i), to_entity, entities);
      if (found)
      {
        entities.push_back(from_entity);
        return true;
      }
    }

    return false;
  }
}

// gets sequence of entities that connect from_entity to to_entity such that
// entity i+1 is on the closure of entity i
// The returned sequence is in order of decreasing dimension, and from_entity
// is the first item in the list.  to_entity does not appear in the sequence,
// the final entity is of one dimension higher.
std::vector<GEPtr> getEntityChain(GEPtr from_entity, GEPtr to_entity)
{
  std::vector<GEPtr> entities;
  bool found = getEntityChain(from_entity, to_entity, entities);
  assertAlways(found, "Unable to find entity chain");

  return entities;
}


Point reclassifyPoint(GEPtr from_entity, const Point& pt_in, GEPtr to_entity)
{
  // allow "reclassifying" on self for convenience
  if (from_entity == to_entity)
    return pt_in;

  assertAlways(from_entity->getDimension() < to_entity->getDimension(), 
          "can only reclassify point from lower dimension entity to higher dimension entity");

  auto entities = getEntityChain(to_entity, from_entity);

  Point pt = pt_in;
  GEPtr prev_entity = from_entity;
  for (auto it = entities.begin(); it != entities.end(); ++it)
  {
    pt = (*it)->computeFromDownCoordinates(prev_entity, pt);
    prev_entity = *it;
  }

  return pt;
}


} // namespace