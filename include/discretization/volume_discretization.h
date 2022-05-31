#ifndef VOLUME_DISCRETIZATION_H
#define VOLUME_DISCRETIZATION_H

#include "ProjectDefs.h"
#include <cassert>
#include "mesh/volume_group.h"
#include "utils/lagrange.h"
#include "utils/quadrature.h"

class VolumeDiscretization
{
  public:
    explicit VolumeDiscretization(const Mesh::VolumeGroup& vol_group,
                                  const Quadrature& quad);
    ArrayType<Real, 4> dxidx;  // num elements x numQuadPtsPerElement x 3 x 3
    ArrayType<Real, 2> detJ; // determinant of dxi/dx, num elements x numQuadPtsPerFace
    ArrayType<Real, 2> detJInv; // determinant of dx/dxi (= 1/detJ)
    const Mesh::VolumeGroup& vol_group;
    Quadrature quad;
    // coordinate to solution interpolation
    LagrangeEvaluatorTPToTP interp_cs_tp_to_tp;
    LagrangeEvaluatorTPFlatToTP interp_cs_flat_to_tp;
    LagrangeEvaluatorTPToTPFlat interp_cs_tp_to_flat;
    LagrangeEvaluatorTPFlatToTPFlat interp_cs_flat_to_flat;

    // coordinate to quadrature interpolation
    Mesh::TensorProductMapper tp_mapper_quad;
    LagrangeEvaluatorTPFlatToTPFlat interp_cq_flat_to_flat;

    // solution to quadrature interpolation
    LagrangeEvaluatorTPFlatToTPFlat interp_sq_flat_to_flat;

    int getIdx() const { return vol_group.getIdx();}

    int getNumElems() const { return vol_group.getNumElems();}

    int getNumSolPtsPerElement() const { return vol_group.getNumSolPtsPerElement();}

    int getNumCoordPtsPerElement() const { return vol_group.getNumCoordPtsPerElement();}

    int getNumQuadPtsPerElement() const { return quad.getNumPoints() * quad.getNumPoints() * quad.getNumPoints(); }

    int_least8_t getElementWeight(int elnum) const { return vol_group.getElementWeight(elnum); }

    // computes coordinates of volume quadrature points in a flat array (num quad points per face x 3)
    template <typename Array2D>
    void getVolumeSolCoords(const Index el, Array2D& quad_coords);

    // computes coordinates of volume quadrature points in a flat array (num quad points per face x 3)
    template <typename Array2D>
    void getVolumeQuadCoords(const Index el, Array2D& quad_coords);
};

using VolDiscPtr = std::shared_ptr<VolumeDiscretization>;


template <typename Array2D>
void VolumeDiscretization::getVolumeSolCoords(const Index el, Array2D& sol_coords)
{
  assert(sol_coords.num_dimensions() == 2);
  assert(sol_coords.shape()[0] == getNumSolPtsPerElement());
  assert(sol_coords.shape()[1] == 3);

  for (int d=0; d < 3; ++d)
  {
    auto coords_d = vol_group.coords[boost::indices[el][range()][d]];
    auto sol_coords_d = sol_coords[boost::indices[range()][d]];
    interp_cs_flat_to_flat.interpolateVals(coords_d, sol_coords_d);
  }
}

template <typename Array2D>
void VolumeDiscretization::getVolumeQuadCoords(const Index el, Array2D& quad_coords)
{
  assert(quad_coords.num_dimensions() == 2);
  assert(quad_coords.shape()[0] == getNumQuadPtsPerElement());
  assert(quad_coords.shape()[1] == 3);

  for (int d=0; d < 3; ++d)
  {
    auto coords_d = vol_group.coords[boost::indices[el][range()][d]];
    auto quad_coords_d = quad_coords[boost::indices[range()][d]];
    interp_cq_flat_to_flat.interpolateVals(coords_d, quad_coords_d);
  }
}

void computeDxidx(const VolumeDiscretization& vol_disc, ArrayType<Real, 4>& dxidx); 

void computeDetJ(const VolumeDiscretization& vol_disc, ArrayType<Real, 2>& detJ, ArrayType<Real, 2>& detJInv);

template <typename Array>
typename Array::element computeDet3x3(const Array& A)
{
  assert(A.num_dimensions() == 2);
  assert(A.shape()[0] == 3);
  assert(A.shape()[1] == 3);

  return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
         A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
         A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}


template <typename ArrayIn, typename ArrayOut>
void computeCofactor(ArrayIn& A, ArrayOut& B)
{

  assert(A.num_dimensions() == 2);
  assert(A.shape()[0] == 3);
  assert(A.shape()[1] == 3);

  assert(B.num_dimensions() == 2);
  assert(B.shape()[0] == 3);
  assert(B.shape()[1] == 3);

  B[0][0] =   A[1][1]*A[2][2] - A[1][2]*A[2][1];
  B[0][1] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
  B[0][2] =   A[1][0]*A[2][0] - A[1][1]*A[2][0];

  B[1][0] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]);
  B[1][1] =   A[0][0]*A[2][2] - A[0][2]*A[2][0];
  B[1][2] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]);

  B[2][0] =   A[0][1]*A[1][2] - A[0][2]*A[1][1];
  B[2][1] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]);
  B[2][2] =   A[0][0]*A[1][1] - A[0][1]*A[1][0];
}

template <typename ArrayIn, typename ArrayOut>
void computeInverse3x3(ArrayIn& A, ArrayOut& B)
{
  auto detA = computeDet3x3(A);
  computeCofactor(A, B);

  // B = (cofactor(A))^T / detA
  B[0][0] /= detA;
  B[1][1] /= detA;
  B[2][2] /= detA;

  auto tmp = B[0][1];
  B[0][1] = B[1][0]/detA;
  B[1][0] = tmp/detA;

  tmp = B[0][2];
  B[0][2] = B[2][0]/detA;
  B[2][0] = tmp/detA;

  tmp = B[1][2];
  B[1][2] = B[2][1]/detA;
  B[2][1] = tmp/detA;
}


// gets xi coordinates of the volume quadrature points in flat array
void getVolumePoints(const VolumeDiscretization& disc, ArrayType<Real, 2>& face_points);


// integrates a scalar quantity over a single element, given the values at the quadrature points
template <typename Array>
typename Array::element integrateVolumeScalar(VolDiscPtr vol, const int el, const Array& vals)
{
  static_assert(std::is_same<typename Array::element, Real>::value, "Array element type must be Real");

  assert(vals.num_dimensions() == 1);
  assert(vals.shape()[0] == static_cast<unsigned int>(vol->getNumQuadPtsPerElement()));

  typename Array::element val = 0;
  int idx = 0;
  for (int i=0; i < vol->quad.getNumPoints(); ++i)
    for (int j=0; j < vol->quad.getNumPoints(); ++j)
      for (int k=0; k < vol->quad.getNumPoints(); ++k)
      {
        auto weight = vol->quad.getWeight(i) * vol->quad.getWeight(j) * vol->quad.getWeight(k);
        val += vals[idx] * weight / vol->detJ[el][idx];
        idx++;
      }
  
  return val;
}



#endif
