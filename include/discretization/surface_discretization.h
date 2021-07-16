#ifndef SURFACE_DISCRETIZATION_H
#define SURFACE_DISCRETIZATION_H

#include "mesh/mesh.h"
#include "utils/quadrature.h"
#include "discretization/volume_discretization.h"
#include <memory>


class SurfaceDiscretization
{
  public:
    SurfaceDiscretization(const Mesh::FaceGroup& face_group, const std::vector<std::shared_ptr<VolumeDiscretization>> volume_discs);

    int getNumFaces() const {return face_group.getNumFaces();}

    int getNumSolPtsPerFace() const { return face_group.getNumSolPtsPerFace();}

    int getNumQuadPtsPerFace() const {return quad.getNumPoints() * quad.getNumPoints();}

    ArrayType<Real, 3> normals;
    const Mesh::FaceGroup& face_group;
    Quadrature quad;
    std::vector<std::shared_ptr<VolumeDiscretization>> volume_discs;
};

using SurfDiscPtr = std::shared_ptr<SurfaceDiscretization>;

void computeNormals(const SurfaceDiscretization& disc, ArrayType<Real, 3>& normals);

void getFacePoints(const SurfaceDiscretization& disc, ArrayType<Real, 3> face_points);


template <typename Array2D, typename ArrayNormalsXi, typename ArrayNormalsX>
void computeNormalNode(const Array2D& dxidx,
                       const ArrayNormalsXi& normals_xi,
                       ArrayNormalsX& normals_x)
{
  for (int i=0; i < 3; ++i)
    normals_x[i] = dxidx[0][i] * normals_xi[0] + 
                   dxidx[1][i] * normals_xi[1] +
                   dxidx[2][i] * normals_xi[2];
}

// integrates a scalar quantity, given the values at the quadrature points
template <typename Array>
typename Array::element integrateFaceScalar(SurfDiscPtr surf, const int face, const Array& vals)
{
  assert(vals.num_dimensions() == 1);
  assert(vals.shape()[0] == static_cast<unsigned int>(surf->getNumQuadPtsPerFace()));

  typename Array::element val = 0;
  int idx = 0;
  for (int i=0; i < surf->quad.getNumPoints(); ++i)
    for (int j=0; j < surf->quad.getNumPoints(); ++j)
    {
      Real area = 0;
      for (int d=0; d < 3; ++d)
        area += surf->normals[face][idx][d] * surf->normals[face][idx][d];
      area = std::sqrt(area);

      val += vals[idx] * area;
    }

  return val;
}

#endif
