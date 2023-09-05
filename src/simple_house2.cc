#include "bounding_box.h"
#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC.h"
#include "discretization/disc_vector.h"
#include "discretization/discretization.h"
#include "discretization/surface_discretization.h"
#include "linear_system/large_matrix_factory.h"
#include "linear_system/large_matrix_petsc.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/air_leakage.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/environment_interface_weather_file.h"
#include "physics/heat/hvac_model.h"
#include "physics/heat/interior_loads.h"
#include "physics/heat/post_processor_environment_data.h"
#include "physics/heat/post_processor_interior.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_position_calculator.h"
#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/source_terms_def.h"
#include "physics/heat/steady_state_temp_calculator.h"
#include "physics/heat/temperature_controller.h"
#include "physics/heat/window_conduction_model.h"
#include "physics/post_processors.h"
#include "physics/post_processor_scheduler.h"
#include "simple_house/simple_house_spec.h"
#include "time_solver/crank_nicolson.h"
#include "time_solver/time_stepper_opts.h"
#include "time_solver/timestep_controller_piecewise.h"
#include "utils/initialization.h"
#include "simple_house/geometry_generator.h"
#include "simple_house/input_parsing.h"

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  for (auto& val : vec)
    os << val << ", ";

  return os;
}

namespace simple_house {

double computeBottomTemp(std::shared_ptr<Heat::EnvironmentInterface> env_interface,
                        std::shared_ptr<Heat::SolarPositionCalculator> solar_position_calc,
                        double t_start, double t_end, double delta_t,
                        const Params& params)
{
  Heat::TarpModel tarp_model(1, 4, 0, params.vertical_vector, params.point_at_zero_altitude, params.met_terrain_index,
                             params.meterological_altitude, params.local_terrain_index);
  Heat::SkyRadiationModel sky_model(params.lawn_emittance, params.vertical_vector);
  Heat::SolarRadiationModel solar_model(params.lawn_absorptivity);

  Heat::SteadyStateTempCaluclator calc(env_interface, solar_position_calc, tarp_model, sky_model, solar_model);
  return calc.calculate(t_start, t_end, delta_t);
}

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
  auto sky_bc = std::make_shared<Heat::SkyRadiationBC>(surf, emittance, vertical_vector);
  auto solar_bc = std::make_shared<Heat::SolarRadiationBC>(surf, absorptivity);

  std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>> bcs = {tarp_bc};
  if (include_solar)
  {
    bcs.push_back(sky_bc);
    bcs.push_back(solar_bc);
  }
  return std::make_shared<Heat::CombinedAirWindSkyNeumannBC>(bcs);
}

std::shared_ptr<DirichletBC> makeBottomBC(SurfDiscPtr surf, Real temp)
{
  if (commRank(MPI_COMM_WORLD) == 0)
    std::cout << "making bottom BC with temp = " << temp << std::endl;
  // for the bottom exterior surface, set temperature = const
  auto f = [=](Real x, Real y, Real z, Real t) { return temp; };
  return makeDirichletBCMMS(surf, f);
}

void createLawnBC(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn, const Params& params)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  int surface_id = generator.getSurfaceId(SurfaceName::Lawn);
  auto surf = disc->getSurfDisc(surface_id);
  Real surface_area = generator.computeLawnSurfaceArea();
  Real perimeter    = generator.computeLawnPerimeter();
                          
  auto bc = createCombinedBC(surf, surface_area, perimeter, params.lawn_roughness_index, 
                             params.lawn_emittance, params.lawn_absorptivity, true);
  heat_eqn->addNeumannBC(bc, true);
  postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorCombinedAirWindSkyBCFlux>("lawn_flux", bc, heat_eqn.get(), MPI_COMM_WORLD));
}

void setExteriorBCs(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn, Real bottom_temp, const Params& params)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  {
    auto surf = disc->getSurfDisc(generator.getSurfaceId(SurfaceName::GroundBottom));
    //auto bc = std::make_shared<Heat::NewtonCooling>(surf, 0.0);  // zero flux BC for the bottom of the
    //                                                        // ground
    //heat_eqn->addNeumannBC(bc, true);
    heat_eqn->addDirichletBC(makeBottomBC(surf, bottom_temp));
  }

  std::vector<SurfaceName> surf_enums = {SurfaceName::SouthExtWall, SurfaceName::EastExtWall, SurfaceName::NorthExtWall,
                                         SurfaceName::WestExtWall, SurfaceName::Roof};
  std::vector<std::string> names = {"south_exterior_wall_flux", "east_exterior_wall_flux", 
                                    "north_exterior_wall_flux", "west_exterior_wall_flux", "roof_flux",
                                    };

  for (int i=0; i <= 4; ++i)
  {
    int surf_id = generator.getSurfaceId(surf_enums[i]);
    auto surf = disc->getSurfDisc(surf_id);
    int direction = generator.getSurfaceDirection(surf_enums[i]);
    Real surface_area = generator.computeExteriorSurfaceArea(direction);
    Real perimeter    = generator.computeExteriorPerimeter(direction);

    //std::cout << "exterior surface " << i << " with name " << names[i] << " has outward normal "
    //          << surf->normals[0][0][0] << ", " << surf->normals[0][0][1] << ", " << surf->normals[0][0][2] << std::endl;

    Real emittance    = i < 4 ? params.ext_wall_emittance    : params.roof_emittance;
    Real absorptivity = i < 4 ? params.ext_wall_absorptivity : params.roof_absorptivity;
    int roughness_index = i < 4 ? params.ext_wall_roughness_index : params.roof_roughness_index;
    auto bc = createCombinedBC(surf, surface_area, perimeter, roughness_index, emittance, absorptivity, true);
    heat_eqn->addNeumannBC(bc, true);
    postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorCombinedAirWindSkyBCFlux>(names[i], bc, heat_eqn.get(), MPI_COMM_WORLD));
  } 

  createLawnBC(generator, heat_eqn, params);
}

void setExteriorWallTempPostProcessors(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  std::vector<SurfaceName> surf_enums = {SurfaceName::GroundBottom, SurfaceName::SouthExtWall, SurfaceName::EastExtWall,
                                         SurfaceName::NorthExtWall, SurfaceName::WestExtWall, SurfaceName::Roof};

  std::vector<std::string> names = {"ground_bottom_temp", "south_exterior_wall_temp", "east_exterior_wall_temp", 
                                    "north_exterior_wall_temp", "west_exterior_wall_temp", "roof_temp"};

  std::vector<SurfDiscPtr> wall_surfaces;
  auto f = [](double val) { return val; };
  for (int i=0; i <= 5; ++i)
  {
    int surf_id = generator.getSurfaceId(surf_enums[i]);
    auto surf = disc->getSurfDisc(surf_id);
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i], f, MPI_COMM_WORLD));

    if (i >= 1 && i <= 4)
      wall_surfaces.push_back(surf);
  }

  postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(wall_surfaces, "exterior_wall_temp", f, MPI_COMM_WORLD));
}

std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>> makeFloorRadiationBCs(GeometryGenerator& generator, DiscPtr disc,
                                                                              const Params& params)
{
  Real floor_area = generator.computeInteriorSurfaceArea(0);
  auto surf = disc->getSurfDisc(generator.getSurfaceId(SurfaceName::Floor));


  std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>> bcs(4);
  std::vector<std::string> names = {"south_rad", "east_rad", "north_rad", "west_rad"};
  for (int i=0; i < 4; ++i)
  {  
    //TODO: get this from geometry generator?
    std::array<Real, 3> window_normal = {0, 0, 0};
    switch(i)
    {
      case 0: {window_normal[1] = -1; break;}
      case 1: {window_normal[0] =  1; break;}
      case 2: {window_normal[1] =  1; break;}
      case 3: {window_normal[0] = -1; break;}
    }

    bcs[i] = std::make_shared<Heat::FloorRadiationBC>(surf, params.window_areas[i], window_normal, params.window_shgc,
                                                      floor_area, params.floor_absorptivity, params.window_shading_angle, names[i]);
  }

  return bcs;
}

void setInteriorBCs(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn, const Params& params,
                    std::vector<NeumannBCPtr>& interior_air_bcs)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  std::vector<SurfaceName> surf_enums = {SurfaceName::Floor, SurfaceName::SouthIntWall, SurfaceName::EastIntWall,
                                          SurfaceName::NorthIntWall, SurfaceName::WestIntWall, SurfaceName::Ceiling};
  std::vector<std::string> names = {"floor_flux", "south_interior_wall_flux", "east_interior_wall_flux", 
                                    "north_interior_wall_flux", "west_interior_wall_flux", "ceiling_flux"};
  for (int i=0; i <= 5; ++i)
  {
    auto surf = disc->getSurfDisc(generator.getSurfaceId(surf_enums[i]));
    int direction = generator.getSurfaceDirection(surf_enums[i]);
    Real surface_area = generator.computeInteriorSurfaceArea(direction);
    Real perimeter    = generator.computeInteriorPerimeter(direction);
    int roughness_index = i == 0 ? params.floor_roughness_index : params.int_wall_roughness_index;

    std::shared_ptr<Heat::AirWindSkyNeumannBC> bc = std::make_shared<Heat::TarpBC>(surf, surface_area, perimeter, 
                                                                                   roughness_index, params.vertical_vector);
    //std::shared_ptr<Heat::AirWindSkyNeumannBC> bc = std::make_shared<Heat::NewtonCoolingFromAir>(surf, 1.0);
    //std::shared_ptr<Heat::AirWindSkyNeumannBC> bc = std::make_shared<Heat::AirWindSkyZeroBC>(surf);  
    interior_air_bcs.push_back(bc);

    if (i == 0)
    {
      auto bcs = makeFloorRadiationBCs(generator, heat_eqn->getDiscretization(), params);
      bcs.push_back(bc);

      auto bc_combined = std::make_shared<Heat::CombinedAirWindSkyNeumannBC>(bcs);
      heat_eqn->addNeumannBC(bc_combined, false);
      postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorCombinedAirWindSkyBCFlux>(names[i], bc_combined, heat_eqn.get(), MPI_COMM_WORLD));
    } else
    {
      heat_eqn->addNeumannBC(bc, false);
      postprocessors->addPostProcessor(std::make_shared<physics::PostProcessorAirWindSkyBCFlux>(names[i], bc, heat_eqn.get(), MPI_COMM_WORLD));
    }
  }  
}

void setInteriorWallTempPostProcessors(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();

  std::vector<SurfaceName> surf_enums = {SurfaceName::Floor, SurfaceName::SouthIntWall, SurfaceName::EastIntWall,
                                          SurfaceName::NorthIntWall, SurfaceName::WestIntWall, SurfaceName::Ceiling};
  std::vector<std::string> names = {"floor_temp", "south_interior_wall_temp", "east_interior_wall_temp", 
                                    "north_interior_wall_temp", "west_interior_wall_temp", "ceiling_temp"};

  std::vector<SurfDiscPtr> wall_surfaces;
  auto f = [](double val) { return val; };
  for (int i=0; i <= 5; ++i)
  {
    auto surf = disc->getSurfDisc(generator.getSurfaceId(surf_enums[i]));
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i], f, MPI_COMM_WORLD));

    if (i >= 1 && i <= 4)
      wall_surfaces.push_back(surf);
  }

  postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(wall_surfaces, "interior_wall_temp", f, MPI_COMM_WORLD));
}

void setUndergroundTempPostProcessors(GeometryGenerator& generator,  std::shared_ptr<Heat::HeatEquationSolar> heat_eqn)
{
  auto disc = heat_eqn->getDiscretization();
  auto postprocessors = heat_eqn->getPostProcessors();
  std::vector<SurfaceName> surf_enums = {SurfaceName::FoundationBottom, SurfaceName::FoundationInsulationBottom, 
                                         SurfaceName::GroundBottomBeneathFoundation};
  std::vector<std::string> names = {"found_bottom_temp", "found_insl_bottom_temp", "ground_bottom_temp_beneath_foundation"};  
  auto f = [](double val) { return val; };
  for (int i=0; i <= 2; ++i)
  {
    auto surf = disc->getSurfDisc(generator.getSurfaceId(surf_enums[i]));
    postprocessors->addPostProcessor(physics::makePostProcessorSurfaceIntegralAverage(surf, names[i], f, MPI_COMM_WORLD));
  }
}

void createSolarThermalSystem(GeometryGenerator& generator, std::shared_ptr<Heat::HeatEquationSolar> heat_eqn, const Params& params)
{
  int vol_group_idx = generator.getVolumeGroupIdx(0, 0, -1);
  auto vol_disc = heat_eqn->getDiscretization()->getVolDisc(vol_group_idx);

  Real delta_z = params.simple_house_spec.foundation_thicknesses[0] / params.simple_house_spec.foundation_numels[0];
  int num_els = std::ceil(params.solar_min_thickness / delta_z);

  double zmin = 0, zmax = 0;
  if (num_els > params.simple_house_spec.foundation_numels[0])
  {
    zmax = 0;
    zmin = -params.simple_house_spec.foundation_thicknesses[0];
  } else
  {
    zmax = -delta_z;
    zmin = -(1 + num_els) * delta_z;    
  }

  utils::BoundingBox box({MIN_REAL, MIN_REAL, zmin}, {MAX_REAL, MAX_REAL, zmax});

  std::shared_ptr<Heat::TemperatureController> controller;
  if (params.solar_thermal_thermostat)
  {
    controller = std::make_shared<Heat::TemperatureControllerHeatQuadratic>(params.interior_air_min_temp, params.interior_air_max_temp);
  } else
  {
    controller = std::make_shared<Heat::TemperatureControllerConstant>();
  }

  auto middle_block_spec = params.simple_house_spec.middle_block;
  Real delta_x = middle_block_spec.xmax - middle_block_spec.xmin;
  Real delta_y = middle_block_spec.ymax - middle_block_spec.ymin;

  Real foundation_volume = delta_x * delta_y * (zmax - zmin)*delta_z;

  auto src_term = std::make_shared<Heat::SourceTermSolarHeating>(
    vol_disc, params.solar_collector_area, params.solar_collector_efficiency, params.solar_collector_emissivity,
    params.solar_collector_normal, foundation_volume, box, controller);

  heat_eqn->addSourceTerm(vol_group_idx, src_term);

  auto postproc = std::make_shared<Heat::SolarThermalPostProcessor>(src_term);
  heat_eqn->getPostProcessors()->addPostProcessor(postproc);
}

/*
timesolvers::TimeStepperOpts getTimeStepperOpts()
{
  timesolvers::TimeStepperOpts opts;
  Real day = 24*60*60;

  //Real delta_t = 300;
  opts.t_start = 0;
  opts.t_end   = 1*day; // 24*60*60;  // 1 day
  //opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerConstant>(delta_t);
  //opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerResidual>(delta_t, 0);

  //std::vector<timesolvers::TimestepControllerPiecewise::TimestepPoint> pts = { {0, 300}, {day, 300}, {2*day, 800},
  //                                                                             //{3*day, 720}, {4*day, 2400}, {5*day, day}, 
  //                                                                             //{6*day, 2*day},
  //                                                                             {2*opts.t_end, 800}
  //                                                                           };

  std::vector<timesolvers::TimestepControllerPiecewise::TimestepPoint> pts = { {0, 300}, {day, 1800},
                                                                               {2*opts.t_end, 1800}
                                                                             };                                                                             

  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerPiecewise>(pts);
  opts.mat_type = linear_system::LargeMatrixType::Petsc;
  opts.nonlinear_abs_tol = 1e-9;
  opts.nonlinear_rel_tol = 1e-8;
  opts.nonlinear_itermax = 30;

  auto matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>();
  matrix_opts->is_structurally_symmetric = true;
  // I'm not sure if the Jacobian is value symmetric because of the
  // nonlinear neumann BCs
  matrix_opts->petsc_opts["ksp_atol"] = "1e-12";
  matrix_opts->petsc_opts["ksp_rtol"] = "1e-50";
  matrix_opts->petsc_opts["ksp_monitor"] = "";
  opts.matrix_opts = matrix_opts;

  opts.solve_auxiliary_equations_combined_system = true;
  opts.precompute_linear_jacobian = true;

  return opts;
}
*/

std::string parseInputFileName(int argc, char* argv[])
{
  if (argc != 2)
    throw std::runtime_error(std::string("incorrect number of argument: usage ") + argv[0] + std::string(" input_file") );

  return argv[1];
}

}

using namespace simple_house;

int main(int argc, char* argv[])
{
  std::string fname = parseInputFileName(argc, argv);

  PetscOptionsSetValue(NULL, "-on_error_abort", "");
  //linear_system::setPetscGlobalOption("log_view", "");


  initialize(argc, argv);

  double t_start_initialize = MPI_Wtime();
  bool am_i_root = commRank(MPI_COMM_WORLD) == 0;

  {
    Params params = parseParams(fname);
    //SimpleHouseSpec spec = createHouseSpec();
    timesolvers::TimeStepperOpts& opts = params.time_stepper_opts; //getTimeStepperOpts();


    GeometryGenerator generator(params.simple_house_spec);
    std::shared_ptr<Mesh::MeshCG> mesh = generator.getMesh();
    if (am_i_root)
    {
      std::cout << "total number of dofs = " << mesh->getNumTotalDofs() << std::endl;
      std::cout << "number of local dofs = " << mesh->getNumDofs() << std::endl;
    }
    mesh->writeVtkFiles("mesh_initial");
    DiscPtr disc = std::make_shared<Discretization>(mesh, 3, 3);


    //Heat::EnvironmentData edata{305, 0, {1, 0, 0}, 250, 750, 0};
    //auto environment_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
    auto environment_interface_variable = std::make_shared<Heat::EnvironmentInterfaceWeatherFile>(params.weather_filename);
    auto environment_interface = environment_interface_variable;
    //auto environment_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(environment_interface_variable->getEnvironmentData(0));

    auto solar_calc_variable = std::make_shared<Heat::SolarPositionCalculatorNaval>(environment_interface_variable->getStartDate(), 7,
                                            Heat::solar::DMSToRadians(35, 6, 24.3576), 
                                            Heat::solar::DMSToRadians(106, 37, 45.0516));

    //auto solar_calc = std::make_shared<Heat::SolarPositionCalculatorConstant>(solar_calc_variable->computePosition(0));
    auto solar_calc = solar_calc_variable;
    
    
    auto air_leakage = std::make_shared<Heat::AirLeakageModelPressure>(params.air_leakage_ach50, 4, generator.computeInteriorVolume(),
                                                                       params.air_cp, params.air_rho);
    auto air_ventilation = std::make_shared<Heat::AirLeakageModelPressure>(0, 4, generator.computeInteriorVolume(), params.air_cp, params.air_rho);
    auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(0);

    Real window_area    = params.window_areas[0] + params.window_areas[1] + params.window_areas[2] + params.window_areas[3];
    auto window_model   = std::make_shared<Heat::WindowConductionModel>(params.window_r_value, window_area);


    Real initial_air_temp = (params.interior_air_min_temp + params.interior_air_max_temp) / 2;
    std::shared_ptr<Heat::HVACModel> hvac_model;
    if (params.disable_hvac)
    {
      hvac_model = std::make_shared<Heat::HVACModelConstant>(0);
    } else 
    {
      if (params.solar_thermal_thermostat)
      {
        hvac_model = std::make_shared<Heat::HVACModelTempOnlyCoolingOnly>(params.interior_air_min_temp, params.interior_air_max_temp,
                                        params.air_rho*params.air_cp, generator.computeInteriorVolume(), params.hvac_restore_time, 2);
      } else
      {
        hvac_model = std::make_shared<Heat::HVACModelTempOnly>(params.interior_air_min_temp, params.interior_air_max_temp,
                                        params.air_rho*params.air_cp, generator.computeInteriorVolume(), params.hvac_restore_time, 2);
      }
    }



    std::vector<NeumannBCPtr> interior_air_bcs;
    auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(params.air_rho * params.air_cp, generator.computeInteriorVolume(),  
                                                    air_leakage, air_ventilation, interior_loads, window_model, hvac_model);

    auto heat_eqn = std::make_shared<Heat::HeatEquationSolar>(disc, solar_calc, environment_interface, air_updator, MPI_COMM_WORLD);
    if (am_i_root)
      std::cout << "initial air temp = " << initial_air_temp << std::endl;
    auto u_aux = makeAuxiliaryEquationsStorage(heat_eqn->getAuxEquations());
    u_aux->getVector(1)[0] = initial_air_temp;

    auto postprocessor_scheduler = std::make_shared<physics::PostProcessorScheduleFixedInterval>(1);
    auto postprocessors = std::make_shared<physics::PostProcessorManager>(postprocessor_scheduler, "simple_house_data.txt", MPI_COMM_WORLD);
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorInterior>(heat_eqn->getAuxEquationsSolar(), air_updator));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorEnvironmentData>(heat_eqn.get()));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorSolarPositionCalculator>(solar_calc));
    heat_eqn->setPostProcessors(postprocessors);

    generator.createVolumeGroups(heat_eqn);

    // make exterior BCS
    double bottom_temp = computeBottomTemp(environment_interface, solar_calc, opts.t_start, opts.t_end, 30, params);
    if (am_i_root)
      std::cout << "bottom temperature = " << bottom_temp << std::endl;
    setExteriorBCs(generator, heat_eqn, bottom_temp, params);
    setExteriorWallTempPostProcessors(generator, heat_eqn);

    // make interior BCS
    setInteriorBCs(generator, heat_eqn, params, interior_air_bcs);
    setInteriorWallTempPostProcessors(generator, heat_eqn);

    setUndergroundTempPostProcessors(generator, heat_eqn);

    createSolarThermalSystem(generator, heat_eqn, params);

    air_updator->setBCs(interior_air_bcs);

    // make postprocessors for air sub model
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorAirLeakage>(air_leakage, "air_leakage"));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorAirLeakage>(air_ventilation, "air_ventilation"));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorWindowConduction>(window_model));
    postprocessors->addPostProcessor(std::make_shared<Heat::PostProcessorInteriorLoads>(interior_loads));


    heat_eqn->initialize();

    // create CN solver

    DiscVectorPtr u = makeDiscVector(disc);
    auto ic_func = [&](const Real& x, const Real& y, const Real& z)
    {
      Real val = 0;
      if (z > 0)
        val = initial_air_temp;
      else
        val = bottom_temp;

      return val;
    };

    u->setFunc(ic_func);
    //u->set(initial_air_temp);  //TODO: maybe set to steady state solution?
    timesolvers::CrankNicolson timesolver(heat_eqn, u, u_aux, opts);

    // run solver
    double t_start_run = MPI_Wtime();
    mesh->getFieldDataManager().attachVector(u, "solution");
    mesh->writeVtkFiles("solution_initial");
    timesolver.solve();
    mesh->writeVtkFiles("solution_final");

    MPI_Barrier(MPI_COMM_WORLD);
    double t_end_run = MPI_Wtime();

    if (am_i_root)
    {
      std::cout << "\n\nFinished simple house run" << std::endl;
      std::cout << "initializing took " << t_start_run - t_start_initialize << " seconds" << std::endl;
      double t_run_elapsed = t_end_run - t_start_run;
      double t_simulated_elapsed = opts.t_end - opts.t_start;
      std::cout << "simulation took " << t_run_elapsed << " seconds to simulate " << t_simulated_elapsed << " seconds"
                << ", which is " << t_simulated_elapsed/t_run_elapsed << "x realtime" << std::endl;
    }

  } // force destructors to run before MPI_Finalize

  finalize();

  return 0;
}