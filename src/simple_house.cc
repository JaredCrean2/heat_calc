#include "discretization/DirichletBC_defs.h"
#include "discretization/disc_vector.h"
#include "discretization/discretization.h"
#include "discretization/surface_discretization.h"
#include "linear_system/large_matrix_factory.h"
#include "linear_system/large_matrix_petsc.h"
#include "mesh/mesh_generator_multi_block.h"
#include "mesh/mesh_input.h"
#include "physics/heat/HeatEquation.h"
#include "mesh/mesh.h"
#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/air_leakage.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/post_processor_environment_data.h"
#include "physics/heat/post_processor_interior.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_position_calculator.h"
#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/window_conduction_model.h"
#include "physics/post_processors.h"
#include "physics/post_processor_scheduler.h"
#include "time_solver/crank_nicolson.h"
#include "utils/initialization.h"

class GeometryGenerator
{
  public:
    GeometryGenerator()
    {
      createMeshSpec();
      m_exterior_lengths = computeExteriorLengths();
      m_generator = std::make_shared<Mesh::MeshGeneratorMultiBlock>(m_spec);
      m_mesh = m_generator->generate();
      collectGeometricFaces();
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
            if (i != 0 || j != 0 || k != 0)
            {
              Heat::VolumeGroupParams params{0, 0, 0};
              if (k == 0)
                params = horizontal_params[std::max(std::abs(i)-1, std::abs(j)-1)];
              else if (k > 0)
                params = ceiling_params[k-1];
              else 
                params = foundation_params[-k - 1];

              heat_eqn->addVolumeGroupParams(params);
            }
    }

    const std::vector<int>& getExteriorGeometricFaces(int face) const { return m_exterior_geometric_faces[face]; }

    const std::vector<int>& getInteriorGeometricFaces(int face) const { return m_interior_geometric_faces[face]; }

    Real computeInteriorVolume()
    {
      auto& spec = m_spec.middle_block;
      Real delta_x = (spec.xmax - spec.xmin);
      Real delta_y = (spec.ymax - spec.ymin);
      Real delta_z = (spec.zmax - spec.zmin);
      return delta_x * delta_y * delta_z;
    }

    Real computeInteriorSurfaceArea(int direction)
    {
      std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                     m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                     m_spec.middle_block.zmax - m_spec.middle_block.zmin};   

      return lengths[direction] * lengths[(direction + 1) % 3];
    }

    Real computeInteriorPerimeterArea(int direction)
    {
      std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                     m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                     m_spec.middle_block.zmax - m_spec.middle_block.zmin};   

      return 2*lengths[direction] + 2*lengths[(direction + 1) % 3];
    }    

    // direction: 0 = xy plane, 1 = yz plane, 2 = xz plane
    Real computeExteriorSurfaceArea(int direction)
    {
      return m_exterior_lengths[direction] * m_exterior_lengths[(direction + 1) % 3];
    }

    Real computeExteriorPerimeter(int direction)
    {
      return 2*m_exterior_lengths[direction] + 2*m_exterior_lengths[(direction + 1) % 3];
    }

  private:
    void createMeshSpec()
    {
      m_spec.create_middle_block = false;
      m_spec.middle_block = Mesh::getMeshSpec(0, 9.23, 0, 15.38, 0, 2.46, 5, 8, 4);
      //m_spec.middle_block = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 5, 8, 4);


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

      m_xrange = {-int(m_spec.numel_minusx.size()), int(m_spec.numel_plusx.size())};
      m_yrange = {-int(m_spec.numel_minusy.size()), int(m_spec.numel_plusy.size())};
      m_zrange = {-int(m_spec.numel_minusz.size()), int(m_spec.numel_plusz.size())};
    }

    std::array<Real, 3> computeExteriorLengths()
    {
      std::array<Real, 3> lengths = {m_spec.middle_block.xmax - m_spec.middle_block.xmin,
                                     m_spec.middle_block.ymax - m_spec.middle_block.ymin,
                                     m_spec.middle_block.zmax - m_spec.middle_block.zmin};
      for (auto& thickness : m_horizontal_thicknesses)
      {
        lengths[0] += thickness;
        lengths[1] += thickness;
      }

      for (auto& thickness : m_ceiling_thicknesses)
        lengths[2] += thickness;

      for (auto& thickness : m_foundation_thicknesses)
        lengths[2] += thickness;   

      return lengths;
    } 

    void collectGeometricFaces()
    {
      m_exterior_geometric_faces.resize(6);
      m_interior_geometric_faces.resize(6);

      for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
        for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
        {
          m_exterior_geometric_faces[0].push_back(m_generator->getSurfaceGeometricId(i, j, m_zrange[0], 0));
          m_exterior_geometric_faces[5].push_back(m_generator->getSurfaceGeometricId(i, j, m_zrange[1], 5));
        }

      for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
        for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
        {
          m_exterior_geometric_faces[1].push_back(m_generator->getSurfaceGeometricId(i, m_yrange[0], k, 1));
          m_exterior_geometric_faces[3].push_back(m_generator->getSurfaceGeometricId(i, m_yrange[1], k, 3));
        }    

      for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
        for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
        {
          m_exterior_geometric_faces[4].push_back(m_generator->getSurfaceGeometricId(m_xrange[0], j, k, 4));
          m_exterior_geometric_faces[2].push_back(m_generator->getSurfaceGeometricId(m_xrange[1], j, k, 2));
        }

      m_interior_geometric_faces[0].push_back(m_generator->getSurfaceGeometricId(0,  0, -1, 5));
      m_interior_geometric_faces[1].push_back(m_generator->getSurfaceGeometricId(0, -1,  0, 3));
      m_interior_geometric_faces[2].push_back(m_generator->getSurfaceGeometricId(1,  0,  0, 4));
      m_interior_geometric_faces[3].push_back(m_generator->getSurfaceGeometricId(0,  1,  0, 1));
      m_interior_geometric_faces[4].push_back(m_generator->getSurfaceGeometricId(-1, 0,  0, 2));
      m_interior_geometric_faces[5].push_back(m_generator->getSurfaceGeometricId(0,  0,  1, 0));
    }

    void createMeshCG()
    {
      std::vector<Mesh::MeshEntityGroupSpec> volume_groups;
      std::vector<Mesh::MeshEntityGroupSpec> bc_groups;
      std::vector<Mesh::MeshEntityGroupSpec> other_groups;

      for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
        for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
          for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
            if (i != 0 || j != 0 || k != 0)
            {
              volume_groups.emplace_back(std::string("volume_group") + std::to_string(i) + std::to_string(j) + std::to_string(k));
              volume_groups.back().addModelEntity(Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i, j, k)));
            }

      for (int i=0; i < 6; ++i)
        bc_groups.emplace_back(std::string("exterior") + std::to_string(i) );

      for (int i=0; i < 6; ++i)
        bc_groups.emplace_back(std::string("interior") + std::to_string(i) );        

      // exterior faces        
      for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
        for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
        {
          bc_groups[0].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i, j, m_zrange[0], 0)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i,  j, m_zrange[0])));
          bc_groups[5].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i, j, m_zrange[1], 5)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i,  j, m_zrange[1])));
        }

      for (int i=m_xrange[0]; i <= m_xrange[1]; ++i)
        for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
        {
          bc_groups[1].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i, m_yrange[0], k, 1)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i, m_yrange[0], k)));
          bc_groups[3].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(i, m_yrange[1], k, 3)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(i, m_yrange[0], k)));
        }    

      for (int j=m_yrange[0]; j <= m_yrange[1]; ++j)
        for (int k=m_zrange[0]; k <= m_zrange[1]; ++k)
        {
          bc_groups[4].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(m_xrange[0], j, k, 4)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(m_xrange[0], j, k)));
          bc_groups[2].addModelEntity(Mesh::ModelEntitySpec(2, m_generator->getSurfaceGeometricId(m_xrange[1], j, k, 2)),
                                      Mesh::ModelEntitySpec(3, m_generator->getVolumeGeometricId(m_xrange[1], j, k)));
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

      m_meshcg = Mesh::createMeshCG(m_mesh, volume_groups, bc_groups, other_groups, 1, 1);
    }

    std::shared_ptr<Mesh::MeshGeneratorMultiBlock> m_generator;
    apf::Mesh2* m_mesh;
    std::array<int, 3> m_xrange;
    std::array<int, 3> m_yrange;
    std::array<int, 3> m_zrange;
    std::array<Real, 3> m_exterior_lengths;
    std::shared_ptr<Mesh::MeshCG> m_meshcg;


    // all the geometric faces on the exterior faces of the geometry,
    // using the standard face numbering
    std::vector<std::vector<int>> m_exterior_geometric_faces;

    std::vector<std::vector<int>> m_interior_geometric_faces;


    Mesh::MultiBlockMeshSpec m_spec;
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
    std::vector<Heat::VolumeGroupParams> foundation_params { {2.25, 2400, 880} };
};

std::shared_ptr<Heat::CombinedAirWindSkyNeumannBC> createCombinedBC(SurfDiscPtr surf, Real surface_area, Real perimeter, int roughness_index,
                                                                    Real emittance, Real absorptivity)
{
  std::array<Real, 3> vertical_vector = {0, 0, 1};
  std::array<Real, 3> point_at_zero_altitude = {0, 0, -1};
  int met_terrain_index = 2;  // rough, wooded country
  Real meterological_altitude = 1836;
  int local_terrain_index = 2;  // rough, wooded country
  auto tarp_bc = std::make_shared<Heat::TarpBC>(surf, surface_area, perimeter, roughness_index, vertical_vector, 
                  point_at_zero_altitude, met_terrain_index, meterological_altitude, local_terrain_index);
  auto sky_bc = std::make_shared<Heat::SkyRadiationBC>(surf, emittance, vertical_vector);
  auto solar_bc = std::make_shared<Heat::SolarRadiationBC>(surf, absorptivity);

  std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>> bcs = {tarp_bc, sky_bc, solar_bc};
  return std::make_shared<Heat::CombinedAirWindSkyNeumannBC>(bcs);
}

std::shared_ptr<DirichletBC> makeBottomBC(SurfDiscPtr surf, Real temp)
{
  // for the bottom exterior surface, set temperature = const
  auto f = [=](Real x, Real y, Real z, Real t) { return temp; };
  return makeDirichletBCMMS(surf, f);
}

void setExteriorBCs(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn, Real bottom_temp)
{
  auto disc = heat_eqn->getDiscretization();
  heat_eqn->addDirichletBC(makeBottomBC(disc->getSurfDisc((0)), bottom_temp)); // TODO: currently using air temperature
  auto postprocessors = heat_eqn->getPostProcessors();
  std::vector<std::string> names = {"east_exterior_wall_flux", "north_exterior_wall_flux", 
                                    "west_exterior_wall_flux", "south_exterior_wall_flux", "roof_flux"};
  for (int i=1; i <= 5; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    int direction = i < 5 ? (i - 1) % 2 : 3;
    Real surface_area = generator.computeExteriorSurfaceArea(direction);
    Real perimeter    = generator.computeExteriorPerimeter(direction);
    // absorptivity and emissivity values for stucco from: https://remdb.nrel.gov/measures.php?gId=12&ctId=216&scId=2374
    // absorptivity and emissivity values for asphalt singles from Medina "Effects of Single Absorptivity, 
    // radient barrier emissivity", International Journal of Energy Research, 2000, 24:665
    Real emittance    = i < 5 ? 0.9 : 0.78;
    Real absorptivity = i < 5 ? 0.75 : 0.78;
    auto bc = createCombinedBC(surf, surface_area, perimeter, 0, emittance, absorptivity);
    heat_eqn->addNeumannBC(bc, true);
    postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorBCFlux>(names[i-1], bc));
  }  
}

void setExteriorWallTempPostProcessors(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  std::vector<std::string> names = {"east_exterior_wall_temp", "north_exterior_wall_temp", 
                                    "west_exterior_wall_temp", "south_exterior_wall_temp", "roof_temp"};
  for (int i=1; i <= 5; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    auto f = [](double val) { return val; };
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i-1], f));
  }  
}


void setInteriorBCs(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  std::array<Real, 3> vertical_vector = {0, 0, 1};
  auto postprocessors = heat_eqn->getPostProcessors();
  std::vector<std::string> names = {"floor_flux", "east_interior_wall_flux", "north_interior_wall_flux", 
                                    "west_interior_wall_flux", "south_interior_wall_flux", "ceiling_flux"};
  for (int i=6; i <= 11; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    int direction = i == 5 || i == 11 ? 3 : (i - 1) % 2;
    Real surface_area = generator.computeExteriorSurfaceArea(direction);
    Real perimeter    = generator.computeExteriorPerimeter(direction);

    //TODO: need to account for solar heating of floor via windows
    auto bc = std::make_shared<Heat::TarpBC>(surf, surface_area, perimeter, 0, vertical_vector);
    heat_eqn->addNeumannBC(bc, false);
    postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorBCFlux>(names[i-6], bc));
  }  
}

void setInteriorWallTempPostProcessors(GeometryGenerator& generator,  std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();
  std::vector<std::string> names = {"floor_average_temp", "east_interior_wall_temp", "north_interior_wall_temp", 
                                    "west_interior_wall_temp", "south_interior_wall_temp", "ceiling_temp"};
  for (int i=6; i <= 11; ++i)
  {
    auto surf = disc->getSurfDisc(i);
    auto f = [](double val) { return val; };
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i-6], f));
  }
}

timesolvers::TimeStepperOpts getTimeStepperOpts()
{
  timesolvers::TimeStepperOpts opts;
  opts.t_start = 0;
  opts.t_end   = 24*60*60;  // 1 day
  opts.delta_t = 60;  // 1 minute
  opts.mat_type = linear_system::LargeMatrixType::Petsc;
  opts.nonlinear_abs_tol = 1e-6;
  opts.nonlinear_rel_tol = 1e-8;
  opts.nonlinear_itermax = 20;

  auto matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>();
  matrix_opts->is_structurally_symmetric = true;
  // I'm not sure if the Jacobian is value symmetric because of the
  // nonlinear neumann BCs
  matrix_opts->petsc_opts["ksp_atol"] = "1e-12";
  matrix_opts->petsc_opts["ksp_rtol"] = "1e-50";
  matrix_opts->petsc_opts["ksp_monitor"] = "";
  opts.matrix_opts = matrix_opts;

  return opts;
}

int main(int argc, char* argv[])
{
  PetscOptionsSetValue(NULL, "-on_error_abort", "");

  initialize(argc, argv);

  {
    GeometryGenerator generator;
    std::shared_ptr<Mesh::MeshCG> mesh = generator.getMesh();
    std::cout << "total number of dofs = " << mesh->getNumTotalDofs() << std::endl;
    std::cout << "number of local dofs = " << mesh->getNumDofs() << std::endl;
    DiscPtr disc = std::make_shared<Discretization>(mesh, 3, 3);


    Heat::SolarPositionCalculator solar_calc(computeJulianDate({1, 1, 2000}), 7,
                                            Heat::solar::DMSToRadians(35, 6, 24.3576), 
                                            Heat::solar::DMSToRadians(106, 37, 45.0516));
    Heat::EnvironmentData edata{305, 0, {1, 0, 0}, 0, 0, 0};
    auto environment_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);

    Real air_rho = 1.007;
    Real air_cp = 1006;
    Real hvac_restore_time = 60 * 5;
    Real air_leakage_ach50 = 7;
    auto air_leakage = std::make_shared<Heat::AirLeakageModelPressure>(air_leakage_ach50, 4, generator.computeInteriorVolume(), air_cp, air_rho);
    auto air_ventilation = std::make_shared<Heat::AirLeakageModelPressure>(0, 4, generator.computeInteriorVolume(), air_cp, air_rho);
    auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(0);
    auto window_model   = std::make_shared<Heat::WindowConductionModel>(1, 0);  //TODO: zero window area

    // air properties from 6000 ft altitude
    Real interior_air_min_temp = 293.15;
    Real interior_air_max_temp = 297.039;
    auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(interior_air_min_temp, interior_air_max_temp,
                                                    air_rho * air_cp, generator.computeInteriorVolume(),  
                                                    air_leakage, air_ventilation, interior_loads, window_model, hvac_restore_time);

    auto heat_eqn = std::make_shared<Heat::HeatEquationSolar>(disc, solar_calc, environment_interface, air_updator);

    auto postprocessor_scheduler = std::make_shared<physics::PostProcessorScheduleFixedInterval>(1);
    auto postprocessors = std::make_shared<physics::PostProcessorManager>(postprocessor_scheduler, "simple_house_data.txt");
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorInterior>(heat_eqn->getAuxEquationsSolar(), air_updator));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorEnvironmentData>(heat_eqn.get()));
    heat_eqn->setPostProcessors(postprocessors);

    generator.createVolumeGroups(heat_eqn);

    // make exterior BCS
    setExteriorBCs(generator, heat_eqn, edata.air_temp);
    setExteriorWallTempPostProcessors(generator, heat_eqn);

    // make interior BCS
    setInteriorBCs(generator, heat_eqn);
    setInteriorWallTempPostProcessors(generator, heat_eqn);

    heat_eqn->initialize();

    // create CN solver

    timesolvers::TimeStepperOpts opts = getTimeStepperOpts();
    DiscVectorPtr u = makeDiscVector(disc);
    u->set(0.5*(interior_air_min_temp + interior_air_max_temp));  //TODO: maybe set to steady state solution?
    timesolvers::CrankNicolson timesolver(heat_eqn, u, opts);

    // run solver
    mesh->getFieldDataManager().attachVector(u, "solution");
    mesh->writeVtkFiles("solution_initial");
    timesolver.solve();
    mesh->writeVtkFiles("solution_final");

  } // force destructors to run before MPI_Finalize

  finalize();

  return 0;
}