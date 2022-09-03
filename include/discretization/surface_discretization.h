#ifndef SURFACE_DISCRETIZATION_H
#define SURFACE_DISCRETIZATION_H

#include "lagrange.h"
#include "mesh/face_group.h"
#include "utils/lagrange2d.h"
#include "utils/quadrature.h"
#include "discretization/volume_discretization.h"
#include <memory>


class SurfaceDiscretization
{
  public:
    SurfaceDiscretization(const Mesh::FaceGroup& face_group,
        const Quadrature& quad,
        const std::vector<VolDiscPtr>& volume_discs);

    int getIdx() const { return face_group.getIdx(); }

    bool getIsDirichlet() const { return face_group.getIsDirichlet(); }

    bool getIsBoundarySurface() const { return face_group.getIsBoundarySurface();}

    int getNumFaces() const { return face_group.getNumFaces();}

    const Mesh::FaceSpec& getFaceSpec(int face) { return face_group.faces[face]; }

    int getNumCoordPtsPerFace() const { return face_group.getNumCoordPtsPerFace(); }

    int getNumSolPtsPerFace() const { return face_group.getNumSolPtsPerFace();}

    int getNumQuadPtsPerFace() const {return quad.getNumPoints() * quad.getNumPoints();}

    int_least8_t getFaceWeight(int face) const { return face_group.getFaceWeight(face); }

    VolDiscPtr getVolDisc(const int face) { return volume_discs[face_group.faces[face].vol_group]; }

    // computes coordinates of face quadrature points in a flat array (num quad points per face x 3)
    template <typename Array2D>
    void getFaceQuadCoords(const Index face, Array2D& quad_coords);

    // gets coordinates of solution nodes for a given face
    template <typename Array2D>
    void getFaceSolCoords(const Index face, Array2D& coords);

    ArrayType<Real, 3> normals;
    const Mesh::FaceGroup& face_group;
    Quadrature quad;
    std::vector<std::shared_ptr<VolumeDiscretization>> volume_discs;

    // volume coords to face quadrature
    std::vector<LagrangeEvaluatorTPToNonTP> interp_vcq_tp;
    std::vector<LagrangeEvaluatorTPFlatToNonTP> interp_vcq_flat;

    // volume coord to face solution
    std::vector<LagrangeEvaluatorTPFlatToNonTP> interp_vcs_flat;

    // volume solution to face quadrature
    std::vector<LagrangeEvaluatorTPFlatToNonTP> interp_vsq_flat;

    ArrayType<LocalIndex, 2> quad_tp_nodemap;
};

using SurfDiscPtr = std::shared_ptr<SurfaceDiscretization>;

template <typename Array2D>
void SurfaceDiscretization::getFaceQuadCoords(const Index face, Array2D& quad_coords)
{
  auto& face_spec = face_group.faces[face];
  auto vol_disc = getVolDisc(face);

  for (int d=0; d < 3; ++d)
  {
    auto coords_vol = vol_disc->vol_group.coords[boost::indices[face_spec.el_group][range()][d]];
    auto coords_face = quad_coords[boost::indices[range()][d]];
    interp_vcq_flat[face_spec.face].interpolateVals(coords_vol, coords_face);    
  }
}


template <typename Array2D>
void SurfaceDiscretization::getFaceSolCoords(const Index face, Array2D& coords)
{
  auto& face_spec = getFaceSpec(face);
  auto& vol_disc = volume_discs[face_spec.vol_group];
  //auto& localidx = vol_disc.nodemap_coord[face][node];
  for (int d=0; d < 3; ++d)
  {
    auto coords_d = vol_disc->vol_group.coords[boost::indices[face_spec.el_group][range()][d]];
    auto coords_out = coords[boost::indices[range()][d]];
    interp_vcs_flat[face_spec.face].interpolateVals(coords_d, coords_out);
  }
}


void computeNormals(const SurfaceDiscretization& disc, ArrayType<Real, 3>& normals);

// gets xi coordinates of surface quadrature points in flat array
void getFaceQuadXi(const SurfaceDiscretization& disc, ArrayType<Real, 3>& face_points);

// get xi coordinate of surface coordinate nodes on flat array
void getFaceSolXi(const SurfaceDiscretization& disc, ArrayType<Real, 3>& face_points);

// get tensor product nodemap associated with getFaceQuadXi
void getFaceTensorProductMap(int npts,
                             ArrayType<LocalIndex, 2>& nodemap);

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

// integrates a scalar quantity over a single face, given the values at the quadrature points
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
      auto weight = surf->quad.getWeight(i) * surf->quad.getWeight(j);

      val += vals[idx++] * weight * area;
    }

  return val;
}

// integrates a vector quantity over a single face, computing v dot n, given the values at the quadrature points
template <typename Array>
typename Array::element integrateFaceVector(SurfDiscPtr surf, const int face, const Array& vals)
{
  using Element = typename Array::element;
  assert(vals.num_dimensions() == 2);
  assert(vals.shape()[0] == static_cast<unsigned int>(surf->getNumQuadPtsPerFace()));
  assert(vals.shape()[1] == 3);

  Element val = 0;
  int idx = 0;
  for (int i=0; i < surf->quad.getNumPoints(); ++i)
    for (int j=0; j < surf->quad.getNumPoints(); ++j)
    {
      Element val_j = 0;
      for (int d=0; d < 3; ++d)
      {
        val_j += vals[idx][d] * surf->normals[face][idx][d];
      }
      auto weight = surf->quad.getWeight(i) * surf->quad.getWeight(j);

      val += weight * val_j;
      idx++;
    }

  return val;
}

template <typename Array>
typename Array::element integrateFaceVector_rev(SurfDiscPtr surf, const int face, Array& vals_bar, Real val_bar)
{
  using Element = typename Array::element;
  assert(vals_bar.num_dimensions() == 2);
  assert(vals_bar.shape()[0] == static_cast<unsigned int>(surf->getNumQuadPtsPerFace()));
  assert(vals_bar.shape()[1] == 3);

  Element val = 0;
  int idx = 0;
  for (int i=0; i < surf->quad.getNumPoints(); ++i)
    for (int j=0; j < surf->quad.getNumPoints(); ++j)
    {
      auto weight = surf->quad.getWeight(i) * surf->quad.getWeight(j);
      Real val_j_bar = weight * val_bar;

      for (int d=0; d < 3; ++d)
      {
        vals_bar[idx][d] += surf->normals[face][idx][d] * val_j_bar;
      }

      idx++;
    }

  return val;
}

#endif
