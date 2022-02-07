#ifndef REFERENCE_ELEMENT_GEOMETRY_TYPES_H
#define REFERENCE_ELEMENT_GEOMETRY_TYPES_H

#include "mesh/reference_element_geometry.h"

namespace reference_element {

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

} // namespace

#endif