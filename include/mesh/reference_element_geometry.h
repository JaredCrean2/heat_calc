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

// The definition of a reference element is in two parts.  The classes
// below define the *local* information ie. face 2 is the face at
// x1 = 1 of the hex.  The struct below defines the *global* information ie.
// the definition of each entity in terms of its downward adjacencies.
// It also includes permutation information that relates the local information
// to the global information.
struct ReferenceElementDef
{
  using Int = int_least8_t;
  Int nverts = 8;
  std::vector<std::array<Int, 2>> edge_verts = { {0, 1}, {1, 2}, {2, 3}, {0, 3},
                                                  {4, 5}, {5, 6}, {6, 7}, {4, 7},
                                                  {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  std::vector<std::vector<Int>> face_edges      = { {0, 1, 2, 3},
                                                    {0, 9, 4, 8},
                                                    {1, 10, 5, 9},
                                                    {2, 11, 6, 10},
                                                    {3, 11, 7, 8},
                                                    {4, 5, 6, 7}};
  std::vector<std::vector<Int>> face_edge_perms = { {0, 0, 0, 1},
                                                    {0, 0, 1, 1},
                                                    {0, 0, 1, 1},
                                                    {0, 0, 1, 1},
                                                    {0, 0, 1, 1},
                                                    {0, 0, 0, 1}};
  std::vector<Int> el_faces      = {0, 1, 2, 3, 4, 5};
  std::vector<Int> el_face_perms = {0, 0, 0, 0, 0, 0};
};


class GeometricEntity
{
  public:
    explicit GeometricEntity(int dim) :
      m_dim(dim)
    {}

    virtual ~GeometricEntity() = default;

    int getDimension() const {return m_dim; }

    int countDown() const { return m_down.size(); }

    GEPtr getDown(int i) { return m_down.at(i); }

    int getDownPerm(int i) { return m_down_perm.at(i); }

    void setDown(GEPtr entity, int perm=0)
    { 
      assertAlways(entity->getDimension() == getDimension() - 1, "Downward adjacency must be one dimension lower");
      m_down.push_back(entity);
      m_down_perm.push_back(perm);
    }

    int countUp() const { return m_up.size(); }

    GEPtr getUp(int i) { return m_up.at(i); }

    void setUp(GEPtr entity) { m_up.push_back(entity); }

    // given a point on down entity i, compute the coordinates of that point
    // in the coordinate frame of the current entity
    Point computeFromDownCoordinates(int i, const Point pt)
    {
      assertAlways(i >= 0 && static_cast<unsigned int>(i) <= m_down.size(), "invalid downward entity");
      validatePoint(pt);
      return computefromDownCoordinatesImpl(i, pt);
    }

    Point computeFromDownCoordinates(GEPtr entity, const Point pt)
    {
      for (size_t i=0; i < m_down.size(); ++i)
        if (entity == m_down[i])
          return computeFromDownCoordinates(i, pt);

      throw std::invalid_argument("entity must be a downward adjacent entity");
    }

  protected:
    virtual Point computefromDownCoordinatesImpl(int i, const Point pt) = 0;

    // checks that pt can exist on an entity of this dimension
    virtual void validatePoint(const Point pt)
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

  private:
    int m_dim;
    std::vector<GEPtr> m_down;
    std::vector<int> m_down_perm;
    std::vector<GEPtr> m_up;
};


class GeometricVertex : public GeometricEntity
{
  public:
    GeometricVertex(int dim) : GeometricEntity(dim) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override
    {
      throw std::runtime_error("GeometricVertex does not have downward entities");
    }
};

// edge is oriented from v0 to v1
class GeometricEdge : public GeometricEntity
{
  public:
    GeometricEdge(int dim) : GeometricEntity(dim) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override
    {
      switch (i)
      {
        case 0: {return {0, 0, 0}; }
        case 1: {return {1, 0, 0}; }
        case 2: {return {1, 1, 0}; }
        case 3: {return {0, 1, 0}; }
        default:
          throw std::invalid_argument("invalid downward index");
      };
    }
};


// edges are oriented counter clockwise, starting from the lower left corner
class GeometricQuad : public GeometricEntity
{
  public:
    GeometricQuad(int dim) : GeometricEntity(dim) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override
    {
      Real xi = getDownPerm(i) == 0 ? pt[0] : 1 - pt[0];
      switch (i)
      {
        case 0: {return {xi, 0, 0}; }
        case 1: {return {1, xi, 0}; }
        case 2: {return {1 - xi, 1, 0}; }
        case 3: {return {0, 1 - xi, 0}; }
        default:
          throw std::invalid_argument("invalid downward index");
      };
    }
};


// edges are oriented counter clockwise, starting from the lower left corner
// Faces are defined by edges, but the verts are given here
//  0: 0, 1, 2, 3 (bottom)
//  1: 0, 1, 5, 4 (y-)
//  2: 1, 2, 6, 5 (x+)
//  3: 2, 3, 7, 6 (y+)
//  4: 0, 3, 7, 4 (x-)
//  5: 4, 5, 6, 7 (top)
class GeometricHex : public GeometricEntity
{
  public:
    GeometricHex(int dim) : GeometricEntity(dim) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override
    {
      // for quads, perm indicates the rotation of the face relative to how
      // the hex expects it to be.  0 = 90 degree, 1 = 180 degrees etc.
      // Thus we have to apply a rotation of -90, -180, etc.
      int perm = getDownPerm(i);
      Real theta = -perm * pi() / 2;

      //TODO: theta only takes on special values
      Real xi1 = std::cos(theta) * pt[0] - std::sin(theta) * pt[1];
      Real xi2 = std::sin(theta) * pt[0] + std::cos(theta) * pt[1];

      //Note: some of the faces have normals that are oriented inwards,
      //      others are oriented outwards
      switch (i)
      {
        case 0: {return {xi1, xi2, 0}; }
        case 1: {return {xi1, 0, xi2}; }
        case 2: {return {1, xi1, xi2}; }
        case 3: {return {1 - xi1, 1, xi2}; }
        case 4: {return {0, xi1, xi2}; }
        case 5: {return {xi1, xi2, 1}; }
        default:
          throw std::invalid_argument("invalid downward index");
      };
    }
};

class ReferenceElement
{
  public:
    explicit ReferenceElement(const ReferenceElementDef& def);

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

// returns true if to_entity is a downward adjacency of from_entity 
// (traversing several dimensions down if necessary)
bool hasDownwardAdjacency(GEPtr from_entity, GEPtr to_entity);

void checkReferenceElementDef(const ReferenceElementDef& def);

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