#include "discretization/volume_discretization.h"
#include "utils/lagrange.h"
#include <iostream>

VolumeDiscretization::VolumeDiscretization(const Mesh::VolumeGroup& _vol_group) :
  vol_group(_vol_group)
{
  computeDxidx(*this, dxidx);
}

void computeDxidx(const VolumeDiscretization& vol_disc, ArrayType<Real, 4>& dxidx)
{
  using range = boost::multi_array_types::index_range;
  dxidx.resize(boost::extents[vol_disc.getNumElems()][vol_disc.getNumSolPtsPerElement()][3][3]);
  const auto& tp_mapper_coord = vol_disc.vol_group.getTPMapperCoord();
  const auto& tp_mapper_sol   = vol_disc.vol_group.getTPMapperSol();
  LagrangeEvaluatorTP lag(tp_mapper_coord.getXi(), tp_mapper_sol.getXi());

  
  ArrayType<Real, 3> coords_tp(tp_mapper_coord.getTPShape());
  ArrayType<Real, 4> dxdxi_tp(tp_mapper_sol.getTPShape()[3]);
  ArrayType<Real, 3> dxdxi_i(boost::extents[vol_disc.getNumSolPtsPerElement()][3][3]);
  auto&& tp_map = tp_mapper_sol.getNodemap();

  for (int i=0; i < vol_disc.getNumElems(); ++i)
  {

    for (int d=0; d < 3; ++d)
    {
      std::cout << "d = " << d << std::endl;

      auto coords_d = vol_disc.vol_group.coords[boost::indices[i][range()][d]];
      tp_mapper_coord.mapToTP(coords_d, coords_tp);
      std::cout << "coords_tp = \n";
      tp_mapper_coord.printTP(std::cout, coords_tp) << std::endl;

      lag.interpolateDerivs(coords_tp, dxdxi_tp);

      for (unsigned int k1=0; k1 < tp_map.shape()[0]; ++k1)
        for (unsigned int k2=0; k2 < tp_map.shape()[1]; ++k2)
          for (unsigned int k3=0; k3 < tp_map.shape()[2]; ++k3)
            for (unsigned int d2=0; d2 < 3; ++d2)
              dxdxi_i[tp_map[k1][k2][k3]][d][d2] = dxdxi_tp[k1][k2][k3][d2];

      std::cout << "finished d loop" << std::endl;
      std::cout << "d = " << d << std::endl;
    }

    for (int j=0; j < vol_disc.getNumSolPtsPerElement(); ++j)
    {
      auto dxdxi_j = dxdxi_i[boost::indices[j][range()][range()]];
      auto dxidx_j = dxidx[boost::indices[i][j][range()][range()]];
      computeInverse3x3(dxdxi_j, dxidx_j);
    }

  }
}

