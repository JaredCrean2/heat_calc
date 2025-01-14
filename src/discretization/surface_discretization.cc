#include "discretization/surface_discretization.h"
#include "mesh/reference_element_geometry.h"
#include "utils/lagrange.h"

SurfaceDiscretization::SurfaceDiscretization(const Mesh::FaceGroup& face_group, const Quadrature& quad, 
    const std::vector<VolDiscPtr>& volume_discs) :
  face_group(face_group),
  quad(quad),
  volume_discs(volume_discs)
{
  this->quad.setDomain(volume_discs[0]->vol_group.ref_el_coord->getXiRange());
  getFaceTensorProductMap(quad.getNumPoints(), quad_tp_nodemap);

  ArrayType<Real, 3> face_points_quad, face_points_sol;
  getFaceQuadXi(*this, face_points_quad);
  getFaceSolXi(*this, face_points_sol);

  auto& tp_mapper_coord = volume_discs[0]->vol_group.getTPMapperCoord();
  auto& tp_mapper_sol   = volume_discs[0]->vol_group.getTPMapperSol();
  for (int i=0; i < face_group.ref_el_coord->getNumFaces(); ++i)
  {
    interp_vcq_tp.emplace_back(tp_mapper_coord.getXi(),
                       face_points_quad[boost::indices[i][range()][range()]]);
    interp_vcq_flat.emplace_back(tp_mapper_coord.getXi(),
                       face_points_quad[boost::indices[i][range()][range()]],
                       tp_mapper_coord.getNodemap());
    interp_vcs_flat.emplace_back(tp_mapper_coord.getXi(),
                       face_points_sol[boost::indices[i][range()][range()]],
                       tp_mapper_coord.getNodemap());
    interp_vsq_flat.emplace_back(tp_mapper_sol.getXi(),
                                 face_points_quad[boost::indices[i][range()][range()]],
                                 tp_mapper_sol.getNodemap());
  }

  computeNormals(*this, normals);

}


void computeNormals(const SurfaceDiscretization& disc, ArrayType<Real, 3>& normals)
{
  int npts = disc.quad.getNumPoints() * disc.quad.getNumPoints();
  normals.resize(boost::extents[disc.getNumFaces()][npts][3]);

  const auto& normals_xi = disc.face_group.ref_el_coord->getNormals();
  const auto& tp_mapper_coord = disc.volume_discs[0]->vol_group.getTPMapperCoord();

  ArrayType<Real, 3> coords_tp(tp_mapper_coord.getTPShape());
  ArrayType<Real, 3> dxdxi_i(boost::extents[npts][3][3]);
  ArrayType<Real, 2> dxidx(boost::extents[3][3]);


  for (int i=0; i < disc.getNumFaces(); ++i)
  {
    auto spec = disc.face_group.faces[i];
    auto normals_i = normals_xi[boost::indices[spec.face][range()]];
    auto vol_disc = disc.volume_discs[spec.vol_group];

    // compute dxdxi
    for (int d=0; d < 3; ++d)
    {
      auto coords_d = vol_disc->vol_group.coords[boost::indices[spec.el_group][range()][d]];
      auto dxdxi_d = dxdxi_i[boost::indices[range()][d][range()]];

      tp_mapper_coord.mapToTP(coords_d, coords_tp);

      disc.interp_vcq_tp[spec.face].interpolateDerivs(coords_tp, dxdxi_d);
    }

    // compute dxidx and normals
    for (int j=0; j < npts; ++j)
    {
      auto dxdxi_j = dxdxi_i[boost::indices[j][range()][range()]];
      auto normals_j = normals[boost::indices[i][j][range()]];
      computeInverse3x3(dxdxi_j, dxidx);
      computeNormalNode(dxidx, normals_i, normals_j);

      // scale by detJ so normal vectors are proprotional to area
      auto detJ = computeDet3x3(dxdxi_j);
      for (int d=0; d < 3; ++d)
        normals_j[d] *= detJ;
    }
  }
}

void getFaceQuadXi(const SurfaceDiscretization& disc, ArrayType<Real, 3>& face_points)
{
  int npts = disc.quad.getNumPoints() * disc.quad.getNumPoints();
  face_points.resize(boost::extents[disc.face_group.ref_el_coord->getNumFaces()][npts][3]);
  std::array<Real, 3> xi_face;

  for (int face=0; face < disc.face_group.ref_el_sol->getNumFaces(); ++face)
  {
    //TODO: this is a potential performance problem
    auto ge_el = disc.face_group.ref_el_sol->getElement();
    auto ge_face = disc.face_group.ref_el_sol->getFace(face);
    for (int i=0; i < disc.quad.getNumPoints(); ++i)
      for (int j=0; j < disc.quad.getNumPoints(); ++j)
      {
        xi_face[0] = disc.quad.getPoint(i);
        xi_face[1] = disc.quad.getPoint(j);
        int idx = disc.quad_tp_nodemap[i][j];
        auto xi_element = reference_element::reclassifyPoint(ge_face, xi_face, ge_el);
        //disc.face_group.ref_el_coord->computeElementXi(face, xi_face, xi_element);
        for (int d=0; d < 3; ++d)
          face_points[face][idx][d] = xi_element[d];
      }
  }
}

void getFaceSolXi(const SurfaceDiscretization& disc, ArrayType<Real, 3>& face_points)
{
  int nfaces    = disc.face_group.ref_el_coord->getNumFaces();
  auto& nodemap = disc.face_group.getFaceNodesSol();
  auto& vol_xi  = disc.face_group.ref_el_sol->getNodeXi();
  face_points.resize(boost::extents[nfaces][disc.getNumSolPtsPerFace()][3]);

  for (int face=0; face < nfaces; ++face)
  {
    for (int j=0; j < disc.getNumSolPtsPerFace(); ++j)
      for (int d=0; d < 3; ++d)
        face_points[face][j][d] = vol_xi[nodemap[face][j]][d];
  }
}


void getFaceTensorProductMap(int npts, ArrayType<LocalIndex, 2>& nodemap)
{
  nodemap.resize(boost::extents[npts][npts]);
  int idx = 0;
  for (int i=0; i < npts; ++i)
    for (int j=0; j < npts; ++j)
      nodemap[i][j] = idx++;
}
