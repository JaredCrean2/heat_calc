#include "discretization/volume_discretization.h"
#include "mesh/mesh.h"
#include "utils/lagrange.h"
#include <iostream>

VolumeDiscretization::VolumeDiscretization(const Mesh::VolumeGroup& vol_group, const Quadrature& _quad) :
  vol_group(vol_group),
  quad(_quad),
  interp_cs_tp_to_tp(vol_group.getTPMapperCoord().getXi(),
                     vol_group.getTPMapperSol().getXi()),
  interp_cs_flat_to_tp(vol_group.getTPMapperCoord().getXi(),
                       vol_group.getTPMapperSol().getXi(),
                       vol_group.getTPMapperCoord().getNodemap()),
  interp_cs_tp_to_flat(vol_group.getTPMapperCoord().getXi(),
                       vol_group.getTPMapperSol().getXi(),
                       vol_group.getTPMapperSol().getNodemap()),
  interp_cs_flat_to_flat(vol_group.getTPMapperCoord().getXi(),
                         vol_group.getTPMapperSol().getXi(),
                         vol_group.getTPMapperCoord().getNodemap(),
                         vol_group.getTPMapperSol().getNodemap())
{
  this->quad.setDomain(vol_group.ref_el_coord->getXiRange());

  tp_mapper_quad = Mesh::TensorProductMapper(quad.getPoints());
  interp_cq_flat_to_flat = LagrangeEvaluatorTPFlatToTPFlat(vol_group.getTPMapperCoord().getXi(),
                                                           tp_mapper_quad.getXi(),
                                                           vol_group.getTPMapperCoord().getNodemap(),
                                                           tp_mapper_quad.getNodemap());
  
  interp_sq_flat_to_flat = LagrangeEvaluatorTPFlatToTPFlat(vol_group.getTPMapperSol().getXi(),
                                                           tp_mapper_quad.getXi(),
                                                           vol_group.getTPMapperSol().getNodemap(),
                                                           tp_mapper_quad.getNodemap());
  computeDxidx(*this, dxidx);
  computeDetJ(*this, detJ);
}

void computeDxidx(const VolumeDiscretization& vol_disc, ArrayType<Real, 4>& dxidx)
{
  using range = boost::multi_array_types::index_range;
  dxidx.resize(boost::extents[vol_disc.getNumElems()][vol_disc.getNumQuadPtsPerElement()][3][3]);

  ArrayType<Real, 3> dxdxi_i(boost::extents[vol_disc.getNumQuadPtsPerElement()][3][3]);

  for (int i=0; i < vol_disc.getNumElems(); ++i)
  {
    for (int d=0; d < 3; ++d)
    {
      auto coords_d = vol_disc.vol_group.coords[boost::indices[i][range()][d]];
      auto dxdxi_d = dxdxi_i[boost::indices[range()][d][range()]];
      vol_disc.interp_cq_flat_to_flat.interpolateDerivs(coords_d, dxdxi_d);
    }

    for (int j=0; j < vol_disc.getNumQuadPtsPerElement(); ++j)
    {
      auto dxdxi_j = dxdxi_i[boost::indices[j][range()][range()]];
      auto dxidx_j = dxidx[boost::indices[i][j][range()][range()]];
      computeInverse3x3(dxdxi_j, dxidx_j);
    }
  }

}


void computeDetJ(const VolumeDiscretization& vol_disc, ArrayType<Real, 2>& detJ)
{
  detJ.resize(boost::extents[vol_disc.getNumElems()][vol_disc.getNumQuadPtsPerElement()]);

  for (int i=0; i < vol_disc.getNumElems(); ++i)
    for (int j=0; j < vol_disc.getNumQuadPtsPerElement(); ++j)
      detJ[i][j] = computeDet3x3(vol_disc.dxidx[boost::indices[i][j][range()][range()]]);
}


void getVolumePoints(const VolumeDiscretization& disc, ArrayType<Real, 2>& vol_points)
{
  int npts = disc.getNumQuadPtsPerElement();
  vol_points.resize(boost::extents[npts][3]);

  int idx = 0;
  for (int i=0; i < disc.quad.getNumPoints(); ++i)
    for (int j=0; j < disc.quad.getNumPoints(); ++j)
      for (int k=0; k < disc.quad.getNumPoints(); ++k)
      {
        vol_points[idx][0] = disc.quad.getPoint(i);
        vol_points[idx][1] = disc.quad.getPoint(j);
        vol_points[idx][2] = disc.quad.getPoint(k);
        idx++;
      }

}