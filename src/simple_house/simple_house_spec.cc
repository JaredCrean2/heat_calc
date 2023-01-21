#include "simple_house/simple_house_spec.h"

namespace simple_house {

Mesh::MultiBlockMeshSpec SimpleHouseSpec::createMeshSpec()
{
  Mesh::MultiBlockMeshSpec spec;

  addBuilding(spec);
  addFoundationInsulation(spec);
  addGround(spec);
  setupMasks(spec);

  return spec;
}

void SimpleHouseSpec::addBuilding(Mesh::MultiBlockMeshSpec& spec)
{
  spec.middle_block = middle_block; // Mesh::getMeshSpec(0, 9.23, 0, 15.38, 0, 2.46, 5, 8, 4);
  //m_spec.middle_block = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 5, 8, 4);

  for (size_t i=0; i < horizontal_thicknesses.size(); ++i)
  {
    spec.numel_plusx.push_back(horizontal_numels[i]);
    spec.numel_plusy.push_back(horizontal_numels[i]);
    spec.numel_minusx.push_back(horizontal_numels[i]);
    spec.numel_minusy.push_back(horizontal_numels[i]);

    spec.thickness_plusx.push_back(horizontal_thicknesses[i]);
    spec.thickness_plusy.push_back(horizontal_thicknesses[i]);
    spec.thickness_minusx.push_back(horizontal_thicknesses[i]);
    spec.thickness_minusy.push_back(horizontal_thicknesses[i]);
  }

  for (size_t i=0; i < ceiling_thicknesses.size(); ++i)
  {
    spec.thickness_plusz.push_back(ceiling_thicknesses[i]);
    spec.numel_plusz.push_back(ceiling_numels[i]);
  }

  for (size_t i=0; i < foundation_thicknesses.size(); ++i)
  {
    spec.thickness_minusz.push_back(foundation_thicknesses[i]);
    spec.numel_minusz.push_back(foundation_numels[i]);
  }
}

void SimpleHouseSpec::addFoundationInsulation(Mesh::MultiBlockMeshSpec& spec)
{
  //TODO: add validation to check lengths of vectors
  for (size_t i=0; i < foundation_insulation_numels.size(); ++i)
  {
    int numel_i = foundation_insulation_numels[i];
    Real thickness_i = foundation_insulation_thicknesses[i];
    spec.numel_plusx.push_back(numel_i);
    spec.numel_plusy.push_back(numel_i);
    spec.numel_minusx.push_back(numel_i);
    spec.numel_minusy.push_back(numel_i);
    spec.numel_minusz.push_back(numel_i);
    
    spec.thickness_plusx.push_back(thickness_i);
    spec.thickness_plusy.push_back(thickness_i);
    spec.thickness_minusx.push_back(thickness_i);
    spec.thickness_minusy.push_back(thickness_i);
    spec.thickness_minusz.push_back(thickness_i);
  }
}

void SimpleHouseSpec::addGround(Mesh::MultiBlockMeshSpec& spec)
{
  if (ground_horizontal_numel > 0)
  {
    spec.numel_plusx.push_back(ground_horizontal_numel);
    spec.numel_plusy.push_back(ground_horizontal_numel);
    spec.numel_minusx.push_back(ground_horizontal_numel);
    spec.numel_minusy.push_back(ground_horizontal_numel);

    spec.thickness_plusx.push_back(ground_horizontal_thickness);
    spec.thickness_plusy.push_back(ground_horizontal_thickness);
    spec.thickness_minusx.push_back(ground_horizontal_thickness);
    spec.thickness_minusy.push_back(ground_horizontal_thickness);    
  }

  if (ground_depth_numel > 0)
  {
    spec.numel_minusz.push_back(ground_depth_numel);
    spec.thickness_minusz.push_back(ground_depth);
  }
}



void SimpleHouseSpec::setupMasks(Mesh::MultiBlockMeshSpec& spec)
{

  int nx = int(spec.numel_plusx.size()) + int(spec.numel_minusx.size()) + 1;
  int ny = int(spec.numel_plusy.size()) + int(spec.numel_minusy.size()) + 1;
  int nz = int(spec.numel_plusz.size()) + int(spec.numel_minusz.size()) + 1;

  int nlayers_ground_thickness = ground_horizontal_numel > 0 ? 1 : 0;
  //int nlayers_ground_depth     = ground_depth_numel > 0 ? 1 : 0;
  int nstart_walls_x = foundation_insulation_numels.size() + nlayers_ground_thickness;
  int nend_walls_x   = nx - foundation_insulation_numels.size() - nlayers_ground_thickness - 1;

  int nstart_walls_y = nstart_walls_x;
  int nend_walls_y   = nend_walls_x;

  int nstart_walls_z = spec.numel_minusz.size();
  ///int nend_walls_z   = nz;

  spec.create_blocks.resize(boost::extents[nx][ny][nz]);
  for (int i=0; i < nx; ++i)
    for (int j=0; j < ny; ++j)
      for (int k=0; k < nz; ++k)
      {
        bool below_ground = k < nstart_walls_z;
        bool val;
        if (below_ground)
        {
          val = true;
        } else
        {
          val = i >= nstart_walls_x && i <= nend_walls_x &&
                j >= nstart_walls_y && j <= nend_walls_y;
        }
        spec.create_blocks[i][j][k] = val;
      }

  spec.create_blocks[spec.numel_minusx.size()][spec.numel_minusy.size()][spec.numel_minusz.size()] = false;
}

}