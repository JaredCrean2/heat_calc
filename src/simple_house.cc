#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC.h"
#include "discretization/disc_vector.h"
#include "discretization/discretization.h"
#include "discretization/surface_discretization.h"
#include "linear_system/large_matrix_factory.h"
#include "linear_system/large_matrix_petsc.h"
#include "mesh/mesh_generator_multi_block.h"
#include "mesh/mesh_input.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/heat/HeatEquation.h"
#include "mesh/mesh.h"
#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/air_leakage.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/environment_interface_weather_file.h"
#include "physics/heat/interior_loads.h"
#include "physics/heat/post_processor_environment_data.h"
#include "physics/heat/post_processor_interior.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_position_calculator.h"
#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/window_conduction_model.h"
#include "physics/post_processors.h"
#include "physics/post_processor_scheduler.h"
#include "simple_house/simple_house_spec.h"
#include "time_solver/crank_nicolson.h"
#include "utils/initialization.h"

namespace simple_house {

// The coordinate system is: North is +y, East is +x, z is into outer space
// The surfaces are 0-5 are the outer surface, 6-11 are the inner surfaces.
// The order is the same as the reference Hex elements: 0 = bottom, 5 = top,
// 1 is xz plane at y = ymin -> S
// 2 is yz plane at x = xmax -> E
// 3 is xz plane at y = ymax -> N
// 4 is yz plane at x = xmin -> W
class GeometryGenerator
{
  public:
    GeometryGenerator(SimpleHouseSpec& spec) :
      m_spec(spec.createMeshSpec())
    {
      createMeshSpec(spec);
      m_exterior_lengths = computeExteriorLengths();
      m_generator = std::make_shared<Mesh::MeshGeneratorMultiBlock>(m_spec);
      m_mesh = m_generator->generate();
      //collectGeometricFaces();
      createMeshCG();
    }

    // creates a MeshCG.  Sets up 1 volume group per mesh block.
    // Also sets up 12 face groups.  The first 6 the exterior surface, the
    // last 6 are the interior surface.  They use the standard ordering:
    // 0 = xy plane at z min
    // 1 = xz plane at y min
    // 2 = yz plane at x max
    // 3 = xz plane at y max
    // 4 = yz plane at x min
    // 5 = xy plane at z max
    std::shared_ptr<Mesh::MeshCG> getMesh() { return m_meshcg; }

    // creates the volume groups and sets the material properties
    void createVolumeGroups(std::shared_ptr<Heat::HeatEquation> heat_eqn)
    {
      for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
        for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
          for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
          {
            std::cout << "for block " << i << ", " << j << ", " << k << ", create_block = " << m_spec.create_blocks[i - m_xrange[0]][j - m_yrange[0]][k - m_zrange[0]] << std::endl;
            if (m_spec.create_blocks[i - m_xrange[0]][j - m_yrange[0]][k - m_zrange[0]])
            {
              Heat::VolumeGroupParams params{0, 0, 0};
              if (k == 0) {
                std::cout << "case 1" << std::endl;
                params = m_house_spec.horizontal_params[std::max(std::abs(i)-1, std::abs(j)-1)];
              } else if (k > 0) {
                std::cout << "case 2" << std::endl;
                params = m_house_spec.ceiling_params[k-1];
              } else if (k <= -1 && k >= -int(m_house_spec.foundation_numels.size()))
              {
                if (i >= (m_xrange[0] + m_num_underground_thicknesses) &&
                    i <= (m_xrange[1] - m_num_underground_thicknesses) &&
                    j >= (m_yrange[0] + m_num_underground_thicknesses) &&
                    j <= (m_yrange[1] + m_num_underground_thicknesses))
                {
                  std::cout << "case 3" << std::endl;
                  params = m_house_spec.foundation_params[-k - 1];
                } else if (i == m_xrange[0] || i == m_xrange[1] || 
                           j == m_yrange[0] || j == m_yrange[1])
                {
                  std::cout << "case 4" << std::endl;
                  params = m_house_spec.ground_params;
                } else  // not in the foundation and not in the ground
                {
                  std::cout << "case 5" << std::endl;
                  params = m_house_spec.foundation_insulation_params;
                }
              } else if (k == m_zrange[0] + 1)
              {
                std::cout << "case 6" << std::endl;
                params = m_house_spec.foundation_insulation_params;
              } else // bottom layer of ground
              {
                std::cout << "case 7" << std::endl;
                params = m_house_spec.ground_params;
              }

              std::cout << "params = " << params.kappa << ", " << params.Cp << ", " << params.rho << std::endl;

              heat_eqn->addVolumeGroupParams(params);
            }
          }
    }

    //const std::vector<int>& getExteriorGeometricFaces(int face) const { return m_exterior_geometric_faces[face]; }

    //const std::vector<int>& getInteriorGeometricFaces(int face) const { return m_interior_geometric_faces[face]; }

    Real computeInteriorVolume()
    {
      auto& spec = m_spec.middle_block;
      Real delta_x = (spec.xmax - spec.xmin);
      Real delta_y = (spec.ymax - spec.ymin);
      Real delta_z = (spec.zmax - spec.zmin);
      return delta_x * delta_y * delta_z;
    }

    int getSurfaceDirection(int surface)
    {
      return m_surface_directions.at(surface);
    }

    Real computeInteriorSurfaceArea(int direction)
    {
      assertAlways(direction >= 0 && direction <= 2, "direction must be in range [0, 2]");
      std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                     m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                     m_spec.middle_block.zmax - m_spec.middle_block.zmin};   

      return lengths[direction] * lengths[(direction + 1) % 3];
    }

    Real computeInteriorPerimeter(int direction)
    {
      assertAlways(direction >= 0 && direction <= 3, "direction must be in range [0, 2]");

      std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                     m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                     m_spec.middle_block.zmax - m_spec.middle_block.zmin};   

      return 2*lengths[direction] + 2*lengths[(direction + 1) % 3];
    }    

    // direction: 0 = xy plane, 1 = yz plane, 2 = xz plane
    Real computeExteriorSurfaceArea(int direction)
    {
      assertAlways(direction >= 0 && direction <= 2, "direction must be in range [0, 2]");
      return m_exterior_lengths[direction] * m_exterior_lengths[(direction + 1) % 3];
    }

    Real computeExteriorPerimeter(int direction)
    {
      assertAlways(direction >= 0 && direction <= 2, "direction must be in range [0, 2]");
      return 2*m_exterior_lengths[direction] + 2*m_exterior_lengths[(direction + 1) % 3];
    }

    Real computeLawnSurfaceArea()
    {
      auto p = computeLawnDimensions();
      return p.first * p.second;
    }

    Real computeLawnPerimeter()
    {
      auto p = computeLawnDimensions();
      return 2*p.first + 2*p.second;
    }

  private:
    void createMeshSpec(SimpleHouseSpec& spec)
    {
      m_house_spec = spec;
      //m_spec = spec.createMeshSpec();
      
      //m_spec.middle_block = Mesh::getMeshSpec(0, 9.23, 0, 15.38, 0, 2.46, 5, 8, 4);
      //m_spec.middle_block = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 5, 8, 4);

/*
      for (size_t i=0; i < m_horizontal_thicknesses.size(); ++i)
      {
        m_spec.numel_plusx.push_back(m_horizontal_numels[i]);
        m_spec.numel_plusy.push_back(m_horizontal_numels[i]);
        m_spec.numel_minusx.push_back(m_horizontal_numels[i]);
        m_spec.numel_minusy.push_back(m_horizontal_numels[i]);

        m_spec.thickness_plusx.push_back(m_horizontal_thicknesses[i]);
        m_spec.thickness_plusy.push_back(m_horizontal_thicknesses[i]);
        m_spec.thickness_minusx.push_back(m_horizontal_thicknesses[i]);
        m_spec.thickness_minusy.push_back(m_horizontal_thicknesses[i]);
      }

      for (size_t i=0; i < m_ceiling_thicknesses.size(); ++i)
      {
        m_spec.thickness_plusz.push_back(m_ceiling_thicknesses[i]);
        m_spec.numel_plusz.push_back(m_ceiling_numels[i]);
      }

      for (size_t i=0; i < m_foundation_thicknesses.size(); ++i)
      {
        m_spec.thickness_minusz.push_back(m_foundation_thicknesses[i]);
        m_spec.numel_minusz.push_back(m_foundation_numels[i]);
      }
*/

      m_xrange = {-int(m_spec.numel_minusx.size()), int(m_spec.numel_plusx.size())};
      m_yrange = {-int(m_spec.numel_minusy.size()), int(m_spec.numel_plusy.size())};
      m_zrange = {-int(m_spec.numel_minusz.size()), int(m_spec.numel_plusz.size())};

      int num_ground_thicknesses = spec.ground_horizontal_numel > 0 ? 1 : 0;
      int num_ground_depths = spec.ground_depth_numel > 0 ? 1 : 0;
      m_num_underground_thicknesses = spec.foundation_insulation_numels.size() + num_ground_thicknesses;
      m_num_underground_depths = spec.foundation_numels.size() + spec.foundation_insulation_numels.size() + num_ground_depths;


/*      
      int nx = m_xrange[1] - m_xrange[0] + 1;
      int ny = m_yrange[1] - m_yrange[0] + 1;
      int nz = m_zrange[1] - m_zrange[0] + 1;
      m_spec.create_blocks.resize(boost::extents[nx][ny][nz]);
      for (int i=0; i < nx; ++i)
        for (int j=0; j < ny; ++j)
          for (int k=0; k < nz; ++k)
            m_spec.create_blocks[i][j][k] = true;

      m_spec.create_blocks[-m_xrange[0]][-m_yrange[0]][-m_zrange[0]] = false;
*/
      //auto spec2 = spec.createMeshSpec();
      //m_spec = spec.createMeshSpec();

      
    }

    std::array<Real, 3> computeExteriorLengths()
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

    //TODO: this appears to be unused
    /*
    void collectGeometricFaces()
    {
      m_exterior_geometric_faces.resize(6);
      m_interior_geometric_faces.resize(6);

      for (int i=m_xrange[0] + m_num_underground_thicknesses; i <= (m_xrange[1] - m_num_underground_thicknesses); ++i)
        for (int j=m_yrange[0] + m_num_underground_thicknesses; j <= (m_yrange[1] - m_num_underground_thicknesses); ++j)
        {
          m_exterior_geometric_faces[0].push_back(m_generator->getSurfaceGeometricId(i, j, m_zrange[0], 0));
          m_exterior_geometric_faces[5].push_back(m_generator->getSurfaceGeometricId(i, j, m_zrange[1], 5));
        }

      for (int i=m_xrange[0] + m_num_underground_thicknesses; i <= (m_xrange[1] - m_num_underground_thicknesses); ++i)
        for (int k=m_zrange[0] + m_num_underground_depths; k <= m_zrange[1]; ++k)
        {
          m_exterior_geometric_faces[1].push_back(m_generator->getSurfaceGeometricId(i, m_yrange[0] + m_num_underground_thicknesses, k, 1));
          m_exterior_geometric_faces[3].push_back(m_generator->getSurfaceGeometricId(i, m_yrange[1] - m_num_underground_thicknesses, k, 3));
        }

      for (int j=m_yrange[0] + m_num_underground_thicknesses; j <= (m_yrange[1] - m_num_underground_thicknesses); ++j)
        for (int k=m_zrange[0] + m_num_underground_depths; k <= m_zrange[1]; ++k)
        {
          m_exterior_geometric_faces[4].push_back(m_generator->getSurfaceGeometricId(m_xrange[0] + m_num_underground_thicknesses, j, k, 4));
          m_exterior_geometric_faces[2].push_back(m_generator->getSurfaceGeometricId(m_xrange[1] - m_num_underground_thicknesses, j, k, 2));
        }

      m_interior_geometric_faces[0].push_back(m_generator->getSurfaceGeometricId(0,  0, -1, 5));
      m_interior_geometric_faces[1].push_back(m_generator->getSurfaceGeometricId(0, -1,  0, 3));
      m_interior_geometric_faces[2].push_back(m_generator->getSurfaceGeometricId(1,  0,  0, 4));
      m_interior_geometric_faces[3].push_back(m_generator->getSurfaceGeometricId(0,  1,  0, 1));
      m_interior_geometric_faces[4].push_back(m_generator->getSurfaceGeometricId(-1, 0,  0, 2));
      m_interior_geometric_faces[5].push_back(m_generator->getSurfaceGeometricId(0,  0,  1, 0));


      for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
        for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
        {
          if (i >= (m_xrange[0] + m_num_underground_thicknesses) && 
              i <= (m_xrange[1] - m_num_underground_thicknesses) &&
              j >= (m_yrange[0] + m_num_underground_thicknesses) &&
              j <= (m_yrange[1] - m_num_underground_thicknesses))
            continue;

          m_lawn_geometric_faces.push_back(m_generator->getSurfaceGeometricId(i, j, -1, 5));
        }


      // the underground surfaces have a zero flux BC, so no need to create a boundary condition at all
    }
    */

    void createMeshCG()
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

      //bc_groups[0].setIsDirichlet(true);

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

    std::pair<Real, Real> computeLawnDimensions()
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

    std::shared_ptr<Mesh::MeshGeneratorMultiBlock> m_generator;
    apf::Mesh2* m_mesh;
    std::array<int, 3> m_xrange;
    std::array<int, 3> m_yrange;
    std::array<int, 3> m_zrange;
    int m_num_underground_thicknesses;
    int m_num_underground_depths;
    std::array<Real, 3> m_exterior_lengths;
    std::shared_ptr<Mesh::MeshCG> m_meshcg;


    // all the geometric faces on the exterior faces of the geometry,
    // using the standard face numbering
    //std::vector<std::vector<int>> m_exterior_geometric_faces;

    //std::vector<std::vector<int>> m_interior_geometric_faces;

    //std::vector<int> m_lawn_geometric_faces;

    std::vector<int> m_surface_directions;

    SimpleHouseSpec m_house_spec;
    Mesh::MultiBlockMeshSpec m_spec;

    /*
    // 3.5 inches of insulation
    std::vector<double> m_horizontal_thicknesses = {0.0897};
    std::vector<int>    m_horizontal_numels      = {5};
    std::vector<Heat::VolumeGroupParams> horizontal_params = { {0.039, 45, 2020} };

    // 6 inches of insulation
    std::vector<double> m_ceiling_thicknesses = {0.1538};
    std::vector<int>    m_ceiling_numels      = {5};
    std::vector<Heat::VolumeGroupParams> ceiling_params = { {0.039, 45, 2020} };


    // 6 inches of concrete
    std::vector<double> m_foundation_thicknesses = {0.1538};
    std::vector<int>    m_foundation_numels      = {5};
    std::vector<Heat::VolumeGroupParams> foundation_params = { {2.25, 2400, 880} };
    */
};

std::shared_ptr<Heat::CombinedAirWindSkyNeumannBC> createCombinedBC(SurfDiscPtr surf, Real surface_area, Real perimeter, int roughness_index,
                                                                    Real emittance, Real absorptivity, bool include_solar)
{
  std::array<Real, 3> vertical_vector = {0, 0, 1};
  std::array<Real, 3> point_at_zero_altitude = {0, 0, -100};
  int met_terrain_index = 2;  // rough, wooded country
  Real meterological_altitude = 1836;
  int local_terrain_index = 2;  // rough, wooded country
  auto tarp_bc = std::make_shared<Heat::TarpBC>(surf, surface_area, perimeter, roughness_index, vertical_vector, 
                  point_at_zero_altitude, met_terrain_index, meterological_altitude, local_terrain_index);
  //auto tarp_bc = std::make_shared<Heat::AirWindSkyZeroBC>(surf);  
  auto sky_bc = std::make_shared<Heat::SkyRadiationBC>(surf, emittance, vertical_vector);
  auto solar_bc = std::make_shared<Heat::SolarRadiationBC>(surf, absorptivity);

  std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>> bcs = {tarp_bc};
  if (include_solar)
  {
    bcs.push_back(sky_bc);
    bcs.push_back(solar_bc);
  }
  return std::make_shared<Heat::CombinedAirWindSkyNeumannBC>(bcs);

  //return std::make_shared<Heat::TarpBC>(surf, surface_area, perimeter, roughness_index, vertical_vector, 
  //                point_at_zero_altitude, met_terrain_index, meterological_altitude, local_terrain_index);
}

std::shared_ptr<DirichletBC> makeBottomBC(SurfDiscPtr surf, Real temp)
{
  // for the bottom exterior surface, set temperature = const
  auto f = [=](Real x, Real y, Real z, Real t) { return temp; };
  return makeDirichletBCMMS(surf, f);
}

void createLawnBC(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  int surface_id = 12;
  auto surf = disc->getSurfDisc(surface_id);
  Real surface_area = generator.computeLawnSurfaceArea();
  Real perimeter    = generator.computeLawnPerimeter();
  //TODO: get real values
  Real emittance    = 1;
  Real absorptivity = 1;
  auto bc = createCombinedBC(surf, surface_area, perimeter, 0, emittance, absorptivity, true);
  heat_eqn->addNeumannBC(bc, true);
  postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorCombinedAirWindSkyBCFlux>("lawn_flux", bc, heat_eqn.get()));
}

void setExteriorBCs(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn, Real bottom_temp)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  //heat_eqn->addDirichletBC(makeBottomBC(disc->getSurfDisc((0)), bottom_temp)); // TODO: currently using air temperature
  {
    auto surf = disc->getSurfDisc(0);
    auto bc = std::make_shared<Heat::NewtonCooling>(surf, 0.0);  // zero flux BC for the bottom of the
                                                            // ground
    //int direction = generator.getSurfaceDirection(0);
    //Real surface_area = generator.computeExteriorSurfaceArea(direction);
    //std::cout << "surface_area = " << surface_area << std::endl;
    //Real perimeter    = generator.computeExteriorPerimeter(direction);
    // absorptivity and emissivity values for stucco from: https://remdb.nrel.gov/measures.php?gId=12&ctId=216&scId=2374
    // absorptivity and emissivity values for asphalt singles from Medina "Effects of Single Absorptivity, 
    // radient barrier emissivity", International Journal of Energy Research, 2000, 24:665
    //Real emittance    = 0.78;
    //Real absorptivity = 0.78;
    //auto bc = createCombinedBC(surf, surface_area, perimeter, 0, emittance, absorptivity, false);
    heat_eqn->addNeumannBC(bc, true);
    //postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorCombinedAirWindSkyBCFlux>("foundation_flux", bc, heat_eqn.get()));    
  }

  std::vector<std::string> names = {"south_exterior_wall_flux", "east_exterior_wall_flux", 
                                    "north_exterior_wall_flux", "west_exterior_wall_flux", "roof_flux",
                                    };

  for (int i=1; i <= 5; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    int direction = generator.getSurfaceDirection(i);
    Real surface_area = generator.computeExteriorSurfaceArea(direction);
    Real perimeter    = generator.computeExteriorPerimeter(direction);
    // absorptivity and emissivity values for stucco from: https://remdb.nrel.gov/measures.php?gId=12&ctId=216&scId=2374
    // absorptivity and emissivity values for asphalt singles from Medina "Effects of Single Absorptivity, 
    // radient barrier emissivity", International Journal of Energy Research, 2000, 24:665
    Real emittance    = i < 5 ? 0.9 : 0.78;
    Real absorptivity = i < 5 ? 0.75 : 0.78;
    auto bc = createCombinedBC(surf, surface_area, perimeter, 0, emittance, absorptivity, true);
    heat_eqn->addNeumannBC(bc, true);
    postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorCombinedAirWindSkyBCFlux>(names[i-1], bc, heat_eqn.get()));
  } 

  createLawnBC(generator, heat_eqn);
}

void setExteriorWallTempPostProcessors(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  std::vector<std::string> names = {"foundation_temp", "south_exterior_wall_temp", "east_exterior_wall_temp", 
                                    "north_exterior_wall_temp", "west_exterior_wall_temp", "roof_temp"};
  auto f = [](double val) { return val; };
  for (int i=0; i <= 5; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i], f));
  }

  std::vector<SurfDiscPtr> wall_surfaces;
  for (int i=1; i <= 4; ++i)
    wall_surfaces.push_back(disc->getSurfDisc(i));

  postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(wall_surfaces, "exterior_wall_temp", f));
}

std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>> makeFloorRadiationBCs(GeometryGenerator& generator, DiscPtr disc, std::array<Real, 4> window_areas)
{
  Real shgc = 0.9;
  Real floor_absorptivity = 0.65; // unfinished concrete
  Real floor_area = generator.computeInteriorSurfaceArea(0);
  auto surf = disc->getSurfDisc(6);


  std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>> bcs(4);
  std::vector<std::string> names = {"south_rad", "east_rad", "north_rad", "west_rad"};
  for (int i=0; i < 4; ++i)
  {
    //FloorRadiationBC(SurfDiscPtr surf, Real window_area, std::array<Real, 3> window_normal,
    //                 Real shgc, Real floor_area, Real floor_absorbtivity) :   
    std::array<Real, 3> window_normal = {0, 0, 0};
    switch(i)
    {
      case 0: {window_normal[1] = -1; break;}
      case 1: {window_normal[0] =  1; break;}
      case 2: {window_normal[1] =  1; break;}
      case 3: {window_normal[0] = -1; break;}
    }

    bcs[i] = std::make_shared<Heat::FloorRadiationBC>(surf, window_areas[i], window_normal, shgc, floor_area, floor_absorptivity, names[i]);
  }

  return bcs;
}

void setInteriorBCs(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn, std::array<Real, 4> window_areas)
{
  auto disc = heat_eqn->getDiscretization();
  std::array<Real, 3> vertical_vector = {0, 0, 1};
  auto postprocessors = heat_eqn->getPostProcessors();
  std::vector<std::string> names = {"floor_flux", "south_interior_wall_flux", "east_interior_wall_flux", 
                                    "north_interior_wall_flux", "west_interior_wall_flux", "ceiling_flux"};
  for (int i=6; i <= 11; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    int direction = generator.getSurfaceDirection(i);
    Real surface_area = generator.computeInteriorSurfaceArea(direction);
    Real perimeter    = generator.computeInteriorPerimeter(direction);

    std::shared_ptr<Heat::AirWindSkyNeumannBC> bc = std::make_shared<Heat::TarpBC>(surf, surface_area, perimeter, 0, vertical_vector);

    if (i == 6)
    {
      auto bcs = makeFloorRadiationBCs(generator, heat_eqn->getDiscretization(), window_areas);
      bcs.push_back(bc);

      auto bc_combined = std::make_shared<Heat::CombinedAirWindSkyNeumannBC>(bcs);
      heat_eqn->addNeumannBC(bc_combined, false);
      postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorCombinedAirWindSkyBCFlux>(names[i-6], bc_combined, heat_eqn.get()));
    } else
    {
      heat_eqn->addNeumannBC(bc, false);
      postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorAirWindSkyBCFlux>(names[i-6], bc, heat_eqn.get()));
    }
  }  
}

void setInteriorWallTempPostProcessors(GeometryGenerator& generator,  std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();
  std::vector<std::string> names = {"floor_temp", "south_interior_wall_temp", "east_interior_wall_temp", 
                                    "north_interior_wall_temp", "west_interior_wall_temp", "ceiling_temp"};
  auto f = [](double val) { return val; };
  for (int i=6; i <= 11; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i-6], f));
  }

  std::vector<SurfDiscPtr> wall_surfaces;
  for (int i=7; i <= 10; ++i)
    wall_surfaces.push_back(disc->getSurfDisc(i));

  postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(wall_surfaces, "interior_wall_temp", f));
}

void setUndergroundTempPostProcessors(GeometryGenerator& generator,  std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();
  std::vector<std::string> names = {"found_bottom_temp", "found_insl_bottom_temp", "ground_bottom_temp"};  
  auto f = [](double val) { return val; };
  for (int i=13; i <= 15; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i-13], f));
  }
}

timesolvers::TimeStepperOpts getTimeStepperOpts()
{
  timesolvers::TimeStepperOpts opts;
  Real delta_t = 450;
  opts.t_start = 0;
  opts.t_end   = 20*365*24*60*60; // 24*60*60;  // 1 day
  //opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerConstant>(delta_t);
  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerResidual>(delta_t, 2);


  opts.mat_type = linear_system::LargeMatrixType::Petsc;
  opts.nonlinear_abs_tol = 1e-9;
  opts.nonlinear_rel_tol = 1e-8;
  opts.nonlinear_itermax = 30;

  auto matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>();
  matrix_opts->is_structurally_symmetric = true;
  // I'm not sure if the Jacobian is value symmetric because of the
  // nonlinear neumann BCs
  matrix_opts->petsc_opts["ksp_atol"] = "1e-14";
  matrix_opts->petsc_opts["ksp_rtol"] = "1e-50";
  matrix_opts->petsc_opts["ksp_monitor"] = "";
  opts.matrix_opts = matrix_opts;

  return opts;
}

SimpleHouseSpec createHouseSpec()
{
  SimpleHouseSpec spec;

  spec.middle_block = Mesh::getMeshSpec(0, 9.23, 0, 15.38, 0, 2.46, 5, 8, 4);

    // 3.5 inches of insulation
  spec.horizontal_thicknesses = {0.0897};
  spec.horizontal_numels      = {5};
  spec.horizontal_params      = { {0.039, 45, 2020} };

  // 6 inches of insulation
  spec.ceiling_thicknesses = {0.1538};
  spec.ceiling_numels      = {5};
  spec.ceiling_params      = { {0.039, 45, 2020} };

  // 6 inches of concrete
  spec.foundation_thicknesses = {0.1538};
  spec.foundation_numels      = {5};
  spec.foundation_params      = { {2.25, 2400, 880} };

  spec.foundation_insulation_numels      = {5};
  spec.foundation_insulation_thicknesses = {0.0897};
  spec.foundation_insulation_params      = {0.039, 45, 2020};

  spec.ground_horizontal_numel     = 10;
  spec.ground_horizontal_thickness = 5;
  
  spec.ground_depth_numel = 10;
  spec.ground_depth       = 10;


  // Values from https://open.library.okstate.edu/rainorshine/chapter/13-2-soil-thermal-properties/
  // Table 13-1, for a 2/3 1/3 mixture of clay and soil organic matter
  spec.ground_params = {2.1, 2200, 1144}; //TODO: get real values

  return spec;
}

}

using namespace simple_house;

int main(int argc, char* argv[])
{
  PetscOptionsSetValue(NULL, "-on_error_abort", "");

  initialize(argc, argv);

  double t_start_initialize = MPI_Wtime();

  {
    SimpleHouseSpec spec = createHouseSpec();

    GeometryGenerator generator(spec);
    std::shared_ptr<Mesh::MeshCG> mesh = generator.getMesh();
    std::cout << "total number of dofs = " << mesh->getNumTotalDofs() << std::endl;
    std::cout << "number of local dofs = " << mesh->getNumDofs() << std::endl;
    DiscPtr disc = std::make_shared<Discretization>(mesh, 3, 3);


    //Heat::EnvironmentData edata{305, 0, {1, 0, 0}, 250, 750, 0};
    //auto environment_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
    auto environment_interface_variable = std::make_shared<Heat::EnvironmentInterfaceWeatherFile>("abq.wea");
    auto environment_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(environment_interface_variable->getEnvironmentData(0));

    auto solar_calc_variable = std::make_shared<Heat::SolarPositionCalculatorNaval>(environment_interface_variable->getJulianDateStart(), 7,
                                            Heat::solar::DMSToRadians(35, 6, 24.3576), 
                                            Heat::solar::DMSToRadians(106, 37, 45.0516));

    auto solar_calc = std::make_shared<Heat::SolarPositionCalculatorConstant>(solar_calc_variable->computePosition(0));
    

    Real air_rho           = 1.007;
    Real air_cp            = 1006;
    Real hvac_restore_time = 60 * 20;
    Real air_leakage_ach50 = 7;
    auto air_leakage = std::make_shared<Heat::AirLeakageModelPressure>(air_leakage_ach50, 4, generator.computeInteriorVolume(), air_cp, air_rho);
    auto air_ventilation = std::make_shared<Heat::AirLeakageModelPressure>(0, 4, generator.computeInteriorVolume(), air_cp, air_rho);
    auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(0);

    std::array<Real, 4> window_areas = {0.557418 * 2, 0.557418 * 2, 0.557418 * 2, 0.557418 * 2};  // 6 sq ft each window, 2 windows per wall
    Real window_area    = window_areas[0] + window_areas[1] + window_areas[2] + window_areas[3];
    Real window_r_value = 3 * 0.1761101838;  // r value converted to SI units
    auto window_model   = std::make_shared<Heat::WindowConductionModel>(window_r_value, window_area);

    // air properties from 6000 ft altitude
    Real interior_air_min_temp = 295; //293.15;
    Real interior_air_max_temp = 295; // 297.039;
    Real initial_air_temp = (interior_air_min_temp + interior_air_max_temp) / 2;
    auto hvac_model = std::make_shared<Heat::HVACModelSwitch>(interior_air_min_temp, interior_air_max_temp, air_rho*air_cp, generator.computeInteriorVolume(), hvac_restore_time);
    //auto hvac_model = std::make_shared<Heat::HVACModelDoubleSpline>(interior_air_min_temp, interior_air_max_temp, air_rho*air_cp, generator.computeInteriorVolume(), hvac_restore_time);


    auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(air_rho * air_cp, generator.computeInteriorVolume(),  
                                                    air_leakage, air_ventilation, interior_loads, window_model, hvac_model);

    auto heat_eqn = std::make_shared<Heat::HeatEquationSolar>(disc, solar_calc, environment_interface, air_updator);
    std::cout << "initial air temp = " << initial_air_temp << std::endl;
    auto u_aux = makeAuxiliaryEquationsStorage(heat_eqn->getAuxEquations());
    u_aux->getVector(1)[0] = initial_air_temp;

    auto postprocessor_scheduler = std::make_shared<physics::PostProcessorScheduleFixedInterval>(1);
    auto postprocessors = std::make_shared<physics::PostProcessorManager>(postprocessor_scheduler, "simple_house_data.txt");
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorInterior>(heat_eqn->getAuxEquationsSolar(), air_updator));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorEnvironmentData>(heat_eqn.get()));
    heat_eqn->setPostProcessors(postprocessors);

    generator.createVolumeGroups(heat_eqn);

    // make exterior BCS
    setExteriorBCs(generator, heat_eqn, environment_interface->getEnvironmentData(0).air_temp );
    setExteriorWallTempPostProcessors(generator, heat_eqn);

    // make interior BCS
    setInteriorBCs(generator, heat_eqn, window_areas);
    setInteriorWallTempPostProcessors(generator, heat_eqn);

    setUndergroundTempPostProcessors(generator, heat_eqn);

    // make postprocessors for air sub model
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorAirLeakage>(air_leakage, "air_leakage"));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorAirLeakage>(air_ventilation, "air_ventilation"));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorWindowConduction>(window_model));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorInteriorLoads>(interior_loads));


    heat_eqn->initialize();

    // create CN solver

    timesolvers::TimeStepperOpts opts = getTimeStepperOpts();
    DiscVectorPtr u = makeDiscVector(disc);
    u->set(initial_air_temp);  //TODO: maybe set to steady state solution?
    timesolvers::CrankNicolson timesolver(heat_eqn, u, u_aux, opts);

    // run solver
    double t_start_run = MPI_Wtime();
    mesh->getFieldDataManager().attachVector(u, "solution");
    mesh->writeVtkFiles("solution_initial");
    timesolver.solve();
    mesh->writeVtkFiles("solution_final");

    double t_end_run = MPI_Wtime();

    std::cout << "\n\nFinished simple house run" << std::endl;
    std::cout << "initializing took " << t_start_run - t_start_initialize << " seconds" << std::endl;
    double t_run_elapsed = t_end_run - t_start_run;
    double t_simulated_elapsed = opts.t_end - opts.t_start;
    std::cout << "simulation took " << t_run_elapsed << " seconds to simulate " << t_simulated_elapsed << " seconds"
              << ", which is " << t_simulated_elapsed/t_run_elapsed << "x realtime" << std::endl;

  } // force destructors to run before MPI_Finalize

  finalize();

  return 0;
}