#ifndef REFERENCE_ELEMENT_GEOMETRY_TYPES_H
#define REFERENCE_ELEMENT_GEOMETRY_TYPES_H

#include "mesh/reference_element_geometry.h"
#include "mesh/reference_element_interface.h"

namespace reference_element {

// The definition of a reference element is in two parts.  The classes
// below define the *local* information ie. face 2 is the face at
// x1 = 1 of the hex.  The struct below defines the *global* information ie.
// the definition of each entity in terms of its downward adjacencies.
// It also includes permutation information that relates the local information
// to the global information.
ReferenceElementDef getStandardReferenceElementDef();


class GeometricVertex : public GeometricEntity
{
  public:
    GeometricVertex(int dim, int id) : GeometricEntity(dim, id) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override;
};

// edge is oriented from v0 to v1
class GeometricEdge : public GeometricEntity
{
  public:
    GeometricEdge(int dim, int id) : GeometricEntity(dim, id) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override;
};


// edges are oriented counter clockwise, starting from the lower left corner
class GeometricQuad : public GeometricEntity
{
  public:
    GeometricQuad(int dim, int id) : GeometricEntity(dim, id) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override;
};



class GeometricHex : public GeometricEntity
{
  public:
    GeometricHex(int dim, int id) : GeometricEntity(dim, id) {}

  protected:
    Point computefromDownCoordinatesImpl(int i, const Point pt) override;
};

} // namespace

#endif