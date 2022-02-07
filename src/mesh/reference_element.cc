#include "mesh/reference_element.h"

namespace Mesh {

// ReferenceElementTensorProduct

ArrayType<Real, 2> ReferenceElementTensorProduct::computeNodeXi() const
{
  auto& nodemap = getTensorProductMap();
  auto& tp_xi   = getTensorProductXi();

  ArrayType<Real, 2> xi(boost::extents[getNumNodes()][3]);
  //static bool is_initialized = false;

  //if (!is_initialized)
 // {
    // TODO: this doesn't work for 2D
    for (int i=0; i < getNumNodesTensorProduct(); ++i)
      for (int j=0; j < getNumNodesTensorProduct(); ++j)
        for (int k=0; k < getNumNodesTensorProduct(); ++k)
        {
            xi[nodemap[i][j][k]][0] = tp_xi[i];
            xi[nodemap[i][j][k]][1] = tp_xi[j];
            xi[nodemap[i][j][k]][2] = tp_xi[k];
        }
  //  is_initialized = true;
  //}

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
  static ArrayType<LocalIndex, 3> degree3(boost::extents[4][4][4]);
  static bool is_initialized = false;

  if (!is_initialized)
  {

    //-------------------------------------------
    // degree 1
    degree1[0][0][0] = 0;
    degree1[1][0][0] = 1;
    degree1[0][1][0] = 3;
    degree1[1][1][0] = 2;
    degree1[0][0][1] = 4;
    degree1[1][0][1] = 5;
    degree1[0][1][1] = 7;
    degree1[1][1][1] = 6;

    //-------------------------------------------
    // degree 2
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


    //-------------------------------------------
    // degree 3
    offset = 8;
    degree3[0][0][0] = 0;
    degree3[3][0][0] = 1;
    degree3[0][3][0] = 3;
    degree3[3][3][0] = 2;
    degree3[0][0][3] = 4;
    degree3[3][0][3] = 5;
    degree3[0][3][3] = 7;
    degree3[3][3][3] = 6;

    // mid-edge nodes (bottom)
    degree3[1][0][0] = 0 + offset;
    degree3[2][0][0] = 1 + offset;
    degree3[3][1][0] = 2 + offset;
    degree3[3][2][0] = 3 + offset;

    degree3[1][3][0] = 4 + offset;
    degree3[2][3][0] = 5 + offset;
    degree3[0][1][0] = 6 + offset;
    degree3[0][2][0] = 7 + offset;

    // mid-edge nodes (vertical)
    degree3[0][0][1] = 8  + offset;
    degree3[0][0][2] = 9  + offset;
    degree3[3][0][1] = 10 + offset;
    degree3[3][0][2] = 11 + offset;

    degree3[3][3][1] = 12 + offset;
    degree3[3][3][2] = 13 + offset;
    degree3[0][3][1] = 14 + offset;
    degree3[0][3][2] = 15 + offset;

    // mid-edge nodes (top)
    degree3[1][0][3] = 16 + offset;
    degree3[2][0][3] = 17 + offset;
    degree3[3][1][3] = 18 + offset;
    degree3[3][2][3] = 19 + offset;

    degree3[1][3][3] = 20 + offset;
    degree3[2][3][3] = 21 + offset;
    degree3[0][1][3] = 22 + offset;
    degree3[0][2][3] = 23 + offset;

    // face nodes
    offset = 8 + 24;
    degree3[1][1][0] = 0 + offset;
    degree3[2][1][0] = 1 + offset;
    degree3[1][2][0] = 2 + offset;
    degree3[2][2][0] = 3 + offset;

    degree3[1][0][1] = 4 + offset;
    degree3[2][0][1] = 5 + offset;
    degree3[1][0][2] = 6 + offset;
    degree3[2][0][2] = 7 + offset;

    degree3[3][1][1] = 8 + offset;
    degree3[3][2][1] = 9 + offset;
    degree3[3][1][2] = 10 + offset;
    degree3[3][2][2] = 11 + offset;

    degree3[1][3][1] = 12 + offset;
    degree3[2][3][1] = 13 + offset;
    degree3[1][3][2] = 14 + offset;
    degree3[2][3][2] = 15 + offset;

    degree3[0][1][1] = 16 + offset;
    degree3[0][2][1] = 17 + offset;
    degree3[0][1][2] = 18 + offset;
    degree3[0][2][2] = 19 + offset;

    degree3[1][1][3] = 20 + offset;
    degree3[2][1][3] = 21 + offset;
    degree3[1][2][3] = 22 + offset;
    degree3[2][2][3] = 23 + offset;

    // interior nodes
    offset = 8 + 24 + 24;
    int idx = 0;
    for (int k=1; k < 3; ++k)
      for (int j=1; j < 3; ++j)
        for (int i=1; i < 3; ++i)
          degree3[i][j][k] = idx++ + offset;

    is_initialized = true;
  }

  switch(degree)
  {
    case 1: return degree1;
    case 2: return degree2;
    case 3: return degree3;
    default:
    {
      auto msg = std::string("unsupported degree: ") + std::to_string(degree);
      throw std::invalid_argument(msg);
    }
  }
}


const std::vector<Real>& HexReferenceElement::getTensorProductXi() const
{
  static std::vector<Real> xi1{0, 1};
  static std::vector<Real> xi2{0, 0.5, 1};
  static std::vector<Real> xi3{0, (-std::sqrt(1.0/5.0) + 1)/2.0, (std::sqrt(1.0/5.0) + 1)/2.0, 1};

  switch(degree)
  {
    case 1: return xi1;
    case 2: return xi2;
    case 3: return xi3;
    default:
    {
      auto msg = std::string("unsupported degree ") + std::to_string(degree);
      throw std::invalid_argument(msg);     
    }
  }
}

const ArrayType<Real, 2>& HexReferenceElement::getNodeXi() const
{
  return m_node_xi;
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
  static HexReferenceElement hexes[4] = { HexReferenceElement(1),  // should implement degree 0
                                          HexReferenceElement(1),
                                          HexReferenceElement(2),
                                          HexReferenceElement(3),
                                        };
  switch(type)
  {
    case apf::Mesh::HEX:
    {
      if (degree < 1 || degree > 3)
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
  // TODO: revise face definitions so they are oriented outwards
  switch (face)
  {
    case 0:
    {
      // oriented v0 - v3 - v2 - v1
      xi_element[0] = xi_face[1];
      xi_element[1] = xi_face[0];
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
      // oriented v1 - v2 - v6 - v5
      xi_element[0] = 1;
      xi_element[1] = xi_face[0];  //TODO: is this backwards?
      xi_element[2] = xi_face[1];
      break;
    }

    case 3:
    {
      // oriented v3 - v7 - v6 - v2
      xi_element[0] = xi_face[1];
      xi_element[1] = 1;
      xi_element[2] = xi_face[0];
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
      xi_element[2] = 1;
      break;
    }

    default:
      throw std::invalid_argument("too many faces");
  }
}


}
