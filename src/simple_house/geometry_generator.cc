#include "simple_house/geometry_generator.h"

namespace simple_house {

// creates the volume groups and sets the material properties
void GeometryGenerator::createVolumeGroups(std::shared_ptr<Heat::HeatEquation> heat_eqn)
{
  int num_horizontal_layers = m_house_spec.horizontal_thicknesses.size();
  int num_foundation_layers = m_house_spec.foundation_thicknesses.size();
  int num_foundation_insulation_layers = m_house_spec.foundation_insulation_thicknesses.size();

  for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
    for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
      for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
      {
        bool inside_house_top_view = i >= -num_horizontal_layers &&
                                     i <=  num_horizontal_layers &&
                                     j >= -num_horizontal_layers &&
                                     j <=  num_horizontal_layers;

        bool inside_foundation_vertical = k <= -1  && k >= -num_foundation_layers   ;     

        bool inside_foundation_insulation_top_view = i >= -(num_horizontal_layers + num_foundation_insulation_layers) &&
                                                     i <=  (num_horizontal_layers + num_foundation_insulation_layers) &&
                                                     j >= -(num_horizontal_layers + num_foundation_insulation_layers) &&
                                                     j <=  (num_horizontal_layers + num_foundation_insulation_layers);

        bool inside_foundation_insulation_side_view = k < 0 && k >= -(num_foundation_layers + num_foundation_insulation_layers);

        bool inside_walls      = k == 0 && inside_house_top_view;
        bool inside_ceiling    = k >= 1;
        bool inside_foundation = inside_foundation_vertical && inside_house_top_view;
        bool inside_foundation_insulation = inside_foundation_insulation_top_view && inside_foundation_insulation_side_view && !inside_foundation;
        bool inside_ground = !inside_walls && !inside_ceiling && !inside_foundation && !inside_foundation_insulation;

        if (m_spec.create_blocks[i - m_xrange[0]][j - m_yrange[0]][k - m_zrange[0]])
        {

          Heat::VolumeGroupParams params{0, 0, 0};
          if (inside_walls)
          {
            params = m_house_spec.horizontal_params[std::max(std::abs(i)-1, std::abs(j)-1)];
          } else if (inside_ceiling)
          {
            params = m_house_spec.ceiling_params[k-1];
          } else if (inside_foundation)
          {
            params = m_house_spec.foundation_params[-k - 1];
          } else if (inside_foundation_insulation)
          {
            if (inside_house_top_view)
            {
              //params = m_house_spec.foundation_insulation_params[std::abs(k) - num_foundation_insulation_layers];
              params = m_house_spec.foundation_insulation_params;
            } else
            {
              //int ilayer = std::abs(i) - num_horizontal_layers;
              //int jlayer = std::abs(j) - num_horizontal_layers;
              //int layer  = std::max(ilayer, jlayer);
              //params = m_house_spec.foundation_insulation_params[layer];
              params = m_house_spec.foundation_insulation_params;
            }
          } else if (inside_ground)
          {
            params = m_house_spec.ground_params;
          } else
            throw std::runtime_error("could not determine block affiliation");

          heat_eqn->addVolumeGroupParams(params);
        }
      }
}


Real GeometryGenerator::computeInteriorVolume()
{
  auto& spec = m_spec.middle_block;
  Real delta_x = (spec.xmax - spec.xmin);
  Real delta_y = (spec.ymax - spec.ymin);
  Real delta_z = (spec.zmax - spec.zmin);
  return delta_x * delta_y * delta_z;
}

int GeometryGenerator::getSurfaceDirection(int surface)
{
  return m_surface_directions.at(surface);
}

Real GeometryGenerator::computeInteriorSurfaceArea(int direction)
{
  assertAlways(direction >= 0 && direction <= 2, "direction must be in range [0, 2]");
  std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                  m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                  m_spec.middle_block.zmax - m_spec.middle_block.zmin};   

  return lengths[direction] * lengths[(direction + 1) % 3];
}

Real GeometryGenerator::computeInteriorPerimeter(int direction)
{
  assertAlways(direction >= 0 && direction <= 3, "direction must be in range [0, 2]");

  std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                  m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                  m_spec.middle_block.zmax - m_spec.middle_block.zmin};   

  return 2*lengths[direction] + 2*lengths[(direction + 1) % 3];
}    

// direction: 0 = xy plane, 1 = yz plane, 2 = xz plane
Real GeometryGenerator::computeExteriorSurfaceArea(int direction)
{
  assertAlways(direction >= 0 && direction <= 2, "direction must be in range [0, 2]");
  return m_exterior_lengths[direction] * m_exterior_lengths[(direction + 1) % 3];
}

Real GeometryGenerator::computeExteriorPerimeter(int direction)
{
  assertAlways(direction >= 0 && direction <= 2, "direction must be in range [0, 2]");
  return 2*m_exterior_lengths[direction] + 2*m_exterior_lengths[(direction + 1) % 3];
}

Real GeometryGenerator::computeLawnSurfaceArea()
{
  auto p = computeLawnDimensions();
  return p.first * p.second;
}

Real GeometryGenerator::computeLawnPerimeter()
{
  auto p = computeLawnDimensions();
  return 2*p.first + 2*p.second;
}

int GeometryGenerator::getSurfaceId(SurfaceName surf_name)
{
  return static_cast<int>(surf_name);
}

void GeometryGenerator::createMeshSpec(SimpleHouseSpec& spec)
{
  m_house_spec = spec;

  m_xrange = {-int(m_spec.numel_minusx.size()), int(m_spec.numel_plusx.size())};
  m_yrange = {-int(m_spec.numel_minusy.size()), int(m_spec.numel_plusy.size())};
  m_zrange = {-int(m_spec.numel_minusz.size()), int(m_spec.numel_plusz.size())};

  int num_ground_thicknesses = spec.ground_horizontal_numel > 0 ? 1 : 0;
  int num_ground_depths = spec.ground_depth_numel > 0 ? 1 : 0;
  m_num_underground_thicknesses = spec.foundation_insulation_numels.size() + num_ground_thicknesses;
  m_num_underground_depths = spec.foundation_numels.size() + spec.foundation_insulation_numels.size() + num_ground_depths;      
}

std::array<Real, 3> GeometryGenerator::computeExteriorLengths()
{
  std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                  m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                  m_spec.middle_block.zmax - m_spec.middle_block.zmin};

  for (int i=0; i < int(m_spec.thickness_plusx.size() - m_num_underground_thicknesses); ++i)
    lengths[0] += m_spec.thickness_plusx[i];

  for (int i=0; i < int(m_spec.thickness_minusx.size() - m_num_underground_thicknesses); ++i)
    lengths[0] += m_spec.thickness_minusx[i];

  for (int i=0; i < int(m_spec.thickness_plusy.size() - m_num_underground_thicknesses); ++i)
    lengths[1] += m_spec.thickness_plusy[i];

  for (int i=0; i < int(m_spec.thickness_minusy.size() - m_num_underground_thicknesses); ++i)
    lengths[1] += m_spec.thickness_minusy[i];

  for (auto t : m_spec.thickness_plusz)
    lengths[2] += t;

  return lengths;
} 


void GeometryGenerator::createMeshCG()
{
  std::vector<Mesh::MeshEntityGroupSpec> volume_groups;
  std::vector<Mesh::MeshEntityGroupSpec> bc_groups;
  std::vector<Mesh::MeshEntityGroupSpec> other_groups;

  for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
    for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
      for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
        if (m_spec.create_blocks[i - m_xrange[0]][j - m_yrange[0]][k - m_zrange[0]])
        {
          volume_groups.emplace_back(std::string("volume_group") + std::to_string(i) + std::to_string(j) + std::to_string(k));
          volume_groups.back().addModelEntity(Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i, j, k)));
        }

  for (int i=0; i < 6; ++i)
    bc_groups.emplace_back(std::string("exterior") + std::to_string(i) );

  for (int i=0; i < 6; ++i)
    bc_groups.emplace_back(std::string("interior") + std::to_string(i) );

  bc_groups.emplace_back(std::string("lawn"));

  // exterior faces        
  for (int i=m_xrange[0] + m_num_underground_thicknesses; i <= (m_xrange[1] - m_num_underground_thicknesses); ++i)
    for (int j=m_yrange[0] + m_num_underground_thicknesses; j <= (m_yrange[1] - m_num_underground_thicknesses); ++j)
    {
      bc_groups[0].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i, j, m_zrange[0], 0)),
                                  Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i,  j, m_zrange[0])));
      bc_groups[5].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i, j, m_zrange[1], 5)),
                                  Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i,  j, m_zrange[1])));
    }

  for (int i=m_xrange[0] + m_num_underground_thicknesses; i <= (m_xrange[1] - m_num_underground_thicknesses); ++i)
    for (int k=m_zrange[0] + m_num_underground_depths; k <= m_zrange[1]; ++k)
    {
      bc_groups[1].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i,
                                                        m_yrange[0] + m_num_underground_thicknesses, k, 1)),
                                  Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i,
                                                        m_yrange[0] + m_num_underground_thicknesses, k)));
      bc_groups[3].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i,
                                                        m_yrange[1] - m_num_underground_thicknesses, k, 3)),
                                  Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i,
                                                        m_yrange[1] - m_num_underground_thicknesses, k)));
    }    

  for (int j=m_yrange[0] + m_num_underground_thicknesses; j <= (m_yrange[1] - m_num_underground_thicknesses); ++j)
    for (int k=m_zrange[0] + m_num_underground_depths; k <= m_zrange[1]; ++k)
    {
      bc_groups[4].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(
                                                              m_xrange[0] + m_num_underground_thicknesses, j, k, 4)),
                                  Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(
                                                            m_xrange[0] + m_num_underground_thicknesses, j, k)));
      bc_groups[2].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(
                                                              m_xrange[1] - m_num_underground_thicknesses, j, k, 2)),
                                  Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(
                                                              m_xrange[1] - m_num_underground_thicknesses, j, k)));
    }   

  bc_groups[0].setIsDirichlet(true);

  bc_groups[0+6].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(0,  0, -1, 5)), 
                                Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(0,  0, -1)));
  bc_groups[1+6].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(0, -1,  0, 3)), 
                                Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(0, -1,  0)));
  bc_groups[2+6].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(1,  0,  0, 4)), 
                                Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(1,  0,  0)));
  bc_groups[3+6].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(0,  1,  0, 1)), 
                                Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(0,  1,  0)));
  bc_groups[4+6].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(-1, 0,  0, 2)), 
                                Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(-1, 0,  0)));
  bc_groups[5+6].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(0,  0,  1, 0)),
                                Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(0,  0,  1)));


  for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
    for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
    {
      if (i >= (m_xrange[0] + m_num_underground_thicknesses) && 
          i <= (m_xrange[1] - m_num_underground_thicknesses) &&
          j >= (m_yrange[0] + m_num_underground_thicknesses) &&
          j <= (m_yrange[1] - m_num_underground_thicknesses))
        continue;

      bc_groups[12].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i, j, -1, 5)),
                                    Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i, j, -1)));
    }


  // bottom of foundation
  other_groups.emplace_back("foundation_bottom");
  other_groups.back().addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(0, 0, -1, 0)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(0, 0, -1)));

  // bottom of foundation insulation
  other_groups.emplace_back("foundation_insulation_bottom");
  other_groups.back().addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(0, 0, -2, 0)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(0, 0, -2)));

  // bottom of foundation
  other_groups.emplace_back("ground_bottom");
  other_groups.back().addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(0, 0, -3, 0)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(0, 0, -3)));     

  m_surface_directions = {0, 1, 2, 1, 2, 0,
                          0, 1, 2, 1, 2, 0};        

  m_meshcg = Mesh::createMeshCG(m_mesh, volume_groups, bc_groups, other_groups, 1, 1);
}

std::pair<Real, Real> GeometryGenerator::computeLawnDimensions()
{
  Real total_thickness_x = m_spec.middle_block.xmax - m_spec.middle_block.xmin;
  Real total_thickness_y = m_spec.middle_block.ymax - m_spec.middle_block.ymin;

  for (auto t : m_spec.thickness_plusx)
    total_thickness_x += t;

  for (auto t : m_spec.thickness_minusx)
    total_thickness_x += t;

  for (auto t : m_spec.thickness_plusy)
    total_thickness_y += t;

  for (auto t : m_spec.thickness_minusy)
    total_thickness_y += t;

  return std::make_pair(total_thickness_x, total_thickness_y);
}        


}