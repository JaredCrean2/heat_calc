#ifndef VOLUME_DISCRETIZATION_H
#define VOLUME_DISCRETIZATION_H

#include "ProjectDefs.h"
#include <cassert>
#include "mesh/mesh.h"
#include "utils/lagrange.h"
#include "utils/quadrature.h"

class VolumeDiscretization
{
  public:
    explicit VolumeDiscretization(const Mesh::VolumeGroup& vol_group,
                                  const Quadrature& quad);
    ArrayType<Real, 4> dxidx;
    const Mesh::VolumeGroup& vol_group;
    Quadrature quad;
    LagrangeEvaluatorTPToTP interp_cs_tp;           // coord to solution
    LagrangeEvaluatorTPFlatToTPFlat interp_cs_flat; // coord to solution

    int getNumElems() const { return vol_group.getNumElems();}

    int getNumSolPtsPerElement() const { return vol_group.getNumSolPtsPerElement();}

    int getNumCoordPtsPerElement() const { return vol_group.getNumCoordPtsPerElement();}

};

using VolDiscPtr = std::shared_ptr<VolumeDiscretization>;

void computeDxidx(const VolumeDiscretization& vol_disc, ArrayType<Real, 4>& dxidx); 

template <typename Array>
typename Array::element computeDet3x3(Array& A)
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


#endif
