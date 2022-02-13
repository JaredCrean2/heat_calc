#ifndef REFERENCE_ELEMENT_GEOMETRY_H
#define REFERENCE_ELEMENT_GEOMETRY_H

#include "ProjectDefs.h"

#include <array>
#include <stdexcept>
#include <cmath>
#include <iostream>

#include "utils/error_handling.h"

namespace reference_element {

inline constexpr double pi() { return std::acos(-1); }

using Point = std::array<Real, 3>;

class GeometricEntity;
using GEPtr = std::shared_ptr<GeometricEntity>;

struct ReferenceElementDef;

class GeometricEntity
{
  public:
    explicit GeometricEntity(int dim, int id) :
      m_dim(dim),
      m_id(id)
    {}

    virtual ~GeometricEntity() = default;

    int getDimension() const {return m_dim; }

    int getId() const { return m_id; }

    int countDown() const { return m_down.size(); }

    GEPtr getDown(int i) { return m_down.at(i); }

    int getDownPerm(int i) { return m_down_perm.at(i); }

    void setDown(GEPtr entity, int perm=0);

    int countUp() const { return m_up.size(); }

    GEPtr getUp(int i) { return m_up.at(i); }

    void setUp(GEPtr entity) { m_up.push_back(entity); }

    // given a point on down entity i, compute the coordinates of that point
    // in the coordinate frame of the current entity
    Point computeFromDownCoordinates(int i, const Point pt);

    Point computeFromDownCoordinates(GEPtr entity, const Point pt);

  protected:
    virtual Point computefromDownCoordinatesImpl(int i, const Point pt) = 0;

    // checks that pt can exist on an entity of this dimension
    virtual void validatePoint(const Point pt);

  private:
    int m_dim;
    int m_id;
    std::vector<GEPtr> m_down;
    std::vector<int> m_down_perm;
    std::vector<GEPtr> m_up;
};


std::vector<GEPtr> getDownward(GEPtr entity, int dim);

// returns true if to_entity is a downward adjacency of from_entity 
// (traversing several dimensions down if necessary)
bool hasDownwardAdjacency(GEPtr from_entity, GEPtr to_entity);

//TODO: this is very slow.  See if there is a way to combine hasDownwardAdjacency
//      into this function, rather than evaluating it so many times
// gets sequence of entities that connect from_entity to to_entity such that
// entity i+1 is on the closure of entity i
// The returned sequence is in order of decreasing dimension, and from_entity
// does is the first item in the list.  to_entity does not appear in the sequence,
// the final entity is of one dimension higher.
std::vector<GEPtr> getEntityChain(GEPtr from_entity, GEPtr to_entity);


Point reclassifyPoint(GEPtr from_entity, const Point& pt_in, GEPtr to_entity);

}  // namespace
#endif