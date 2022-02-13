#include "mesh/reference_element_geometry_types.h"

namespace reference_element {

ReferenceElementDef getStandardReferenceElementDef()
{
  using Int = ReferenceElementDef::Int;

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

  std::vector<Int> verts_to_apf = {3, 0, 1, 2, 7, 4, 5, 6};
  std::vector<Int> edges_to_apf = {3, 0, 1, 2, 11, 8, 9, 10, 7, 4, 6, 7};
  std::vector<Int> faces_to_apf = {0, 4, 1, 2, 3, 5};

  ReferenceElementDef def{nverts, edge_verts, face_edges, face_edge_perms, el_faces, el_face_perms,
          verts_to_apf, edges_to_apf, faces_to_apf}; 

  return def;    
};


Point GeometricVertex::computefromDownCoordinatesImpl(int i, const Point pt)
{
  throw std::runtime_error("GeometricVertex does not have downward entities");
}


Point GeometricEdge::computefromDownCoordinatesImpl(int i, const Point pt)
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


Point GeometricQuad::computefromDownCoordinatesImpl(int i, const Point pt)
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

// edges are oriented counter clockwise, starting from the lower left corner
// Faces are defined by edges, but the verts are given here
//  0: 0, 1, 2, 3 (bottom)
//  1: 0, 1, 5, 4 (y-)
//  2: 1, 2, 6, 5 (x+)
//  3: 2, 3, 7, 6 (y+)
//  4: 0, 3, 7, 4 (x-)
//  5: 4, 5, 6, 7 (top)
Point GeometricHex::computefromDownCoordinatesImpl(int i, const Point pt)
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


} // namespace