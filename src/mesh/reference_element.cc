#include "mesh/reference_element.h"

namespace Mesh {

// ReferenceElementTensorProduct

const ArrayType<Real, 2>& ReferenceElementTensorProduct::getNodeXi() const
{
  auto& nodemap = getTensorProductMap();
  auto& tp_xi   = getTensorProductXi();
  static ArrayType<Real, 2> xi(boost::extents[getNumNodes()][3]);
  static bool is_initialized = false;

  if (!is_initialized)
  {
    // TODO: this doesn't work for 2D
    for (int i=0; i < getNumNodesTensorProduct(); ++i)
      for (int j=0; j < getNumNodesTensorProduct(); ++j)
        for (int k=0; k < getNumNodesTensorProduct(); ++k)
        {
            xi[nodemap[i][j][k]][0] = tp_xi[i];
            xi[nodemap[i][j][k]][1] = tp_xi[j];
            xi[nodemap[i][j][k]][2] = tp_xi[k];
        }
    is_initialized = true;
  }

  return xi;
}

int ReferenceElementTensorProduct::getNumNodes() const
{
  const auto& shape = getTensorProductMap().shape();
  return shape[0] * shape[1] * shape[2];
}


// HexReferenceElement
const ArrayType<LocalIndex, 3>& HexReferenceElement::getTensorProductMap() const
{
  static ArrayType<LocalIndex, 3> degree1(boost::extents[2][2][2]);
  static ArrayType<LocalIndex, 3> degree2(boost::extents[3][3][3]);
  static bool is_initialized = false;

  if (!is_initialized)
  {
    degree1[0][0][0] = 0;
    degree1[1][0][0] = 1;
    degree1[0][1][0] = 3;
    degree1[1][1][0] = 2;
    degree1[0][0][1] = 4;
    degree1[1][0][1] = 5;
    degree1[0][1][1] = 7;
    degree1[1][1][1] = 6;

    LocalIndex offset = 8;
    degree2[0][0][0] = 0;
    degree2[2][0][0] = 1;
    degree2[0][2][0] = 3;
    degree2[2][2][0] = 2;
    degree2[0][0][2] = 4;
    degree2[2][0][2] = 5;
    degree2[0][2][2] = 7;
    degree2[2][2][2] = 6;

    // mid-edge nodes (bottom)
    degree2[1][0][0] = 0 + offset;
    degree2[2][1][0] = 1 + offset;
    degree2[1][2][0] = 2 + offset;
    degree2[0][1][0] = 3 + offset;
    // mid-edge nodes (vertical)
    degree2[0][0][1] = 4 + offset;
    degree2[2][0][1] = 5 + offset;
    degree2[2][2][1] = 6 + offset;
    degree2[0][2][1] = 7 + offset;
    // mid-edge nodes (top)
    degree2[1][0][2] = 8  + offset;
    degree2[2][1][2] = 9  + offset;
    degree2[1][2][2] = 10 + offset;
    degree2[0][1][2] = 11 + offset;
    // face nodes
    offset = 8 + 12;
    degree2[1][1][0] = 0 + offset;
    degree2[2][1][1] = 2 + offset;
    degree2[1][2][1] = 3 + offset;
    degree2[0][1][1] = 4 + offset;
    degree2[1][0][1] = 1 + offset;
    degree2[1][1][2] = 5 + offset;
    // interior nodes
    offset = 8 + 12 + 6;
    degree2[1][1][1] = 0 + offset;

    is_initialized = true;
  }

  if (degree == 1)
    return degree1;
  else if (degree == 2)
    return degree2;
  else
  {
    auto msg = std::string("unsupported degree: ") + std::to_string(degree);
    throw std::invalid_argument(msg);
  }
}


const std::vector<Real>& HexReferenceElement::getTensorProductXi() const
{
  static std::vector<Real> xi1{0, 1};
  static std::vector<Real> xi2{0, 0.5, 1};

  if (degree == 1)
    return xi1;
  else if (degree == 2)
    return xi2;
  else
  {
    auto msg = std::string("unsupported degree ") + std::to_string(degree);
    throw std::invalid_argument(msg);
  }
}


const ArrayType<Real, 2>& HexReferenceElement::getNormals() const
{
  // convention: v3 to v1 is xi1, v3 to v2 is xi2, v3 to v7 is xi3
  // convention: v0 to v1 to xi1, v0 to v3 is xi2, v0 to v4 is xi3
  static ArrayType<Real, 2> normals_xi(boost::extents[6][3]);
  static bool is_initialized = false;

  if (!is_initialized)
  {
    normals_xi[0][0] =  0; normals_xi[0][1] =  0; normals_xi[0][2] = -1;
    normals_xi[1][0] =  0; normals_xi[1][1] = -1; normals_xi[1][2] =  0;
    normals_xi[2][0] =  1; normals_xi[2][1] =  0; normals_xi[2][2] =  0;
    normals_xi[3][0] =  0; normals_xi[3][1] =  1; normals_xi[3][2] =  0;
    normals_xi[4][0] = -1; normals_xi[4][1] =  0; normals_xi[4][2] =  0;
    normals_xi[5][0] =  0; normals_xi[5][1] =  0; normals_xi[5][2] =  1;
  }

  return normals_xi;
}


ReferenceElement* getReferenceElement(apf::Mesh::Type type, const int degree)
{
  static HexReferenceElement hexes[3] = { HexReferenceElement(1),  // should implement degree 0
                                          HexReferenceElement(1),
                                          HexReferenceElement(2)
                                        };
  switch(type)
  {
    case apf::Mesh::HEX:
    {
      if (degree < 1 || degree > 2)
        throw std::invalid_argument("unsupported degree");

      return  &(hexes[degree]);
    }

    default: 
      throw std::invalid_argument("unsupported element");
  }
}


void HexReferenceElement::computeElementXi(const int face, const Real* xi_face, Real* xi_element) const
{
  // for all faces, xi1 is from 1st to 2nd vertex, 1st to 4th vertex
  switch (face)
  {
    case 0:
    {
      // oriented v0 - v1 - v2 - v3
      xi_element[0] = xi_face[0];
      xi_element[1] = xi_face[1];
      xi_element[2] = 0;
      break;
    }

    case 1:
    {
      // oriented v0 - v1 - v5 - v4
      xi_element[0] = xi_face[0];
      xi_element[1] = 0;
      xi_element[2] = xi_face[1];
      break;
    }

    case 2:
    {
      // oriented v1 - v5 - v6 - v2
      xi_element[0] = 1;
      xi_element[1] = xi_face[2];
      xi_element[2] = xi_face[1];
      break;
    }

    case 3:
    {
      // oriented v3 - v2 - v6 - v7
      xi_element[0] = xi_face[0];
      xi_element[1] = 1;
      xi_element[2] = xi_face[1];
      break;
    }

    case 4:
    {
      // oriented v0 - v4 - v7 - v3
      xi_element[0] = 0;
      xi_element[1] = xi_face[1];
      xi_element[2] = xi_face[0];
      break;
    }

    case 5:
    {
      // oriented v4 - v5 - v6 - v7
      xi_element[0] = xi_face[0];
      xi_element[1] = xi_face[1];
      xi_element[3] = 1;
      break;
    }

    default:
      throw std::invalid_argument("too many faces");
  }
}


}
