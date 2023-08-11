#include "simple_house/input_parsing.h"

namespace simple_house {

const std::map<std::string, std::string>& getInputFileDefaults()
{
  const static std::map<std::string, std::string> input_file_defaults =
  {
    std::make_pair("met_terrain_index", "2"),
    std::make_pair("meterological_altitude", "1836"),
    std::make_pair("local_terrain_index", "2"),

    std::make_pair("ext_wall_roughness_index", "0"),
    std::make_pair("ext_wall_emittance", "0.9"),
    std::make_pair("ext_wall_absorptivity", "0.75"),

    std::make_pair("roof_roughness_index", "0"),
    std::make_pair("roof_emittance", "0.78"),
    std::make_pair("roof_absorptivity", "0.78"),

    std::make_pair("lawn_roughness_index", "0"),
    std::make_pair("lawn_emittance", "0.95"),
    std::make_pair("lawn_absorptivity", "0.8"), 

    std::make_pair("int_wall_roughness_index", "0"),

    std::make_pair("air_rho", "1.007"),
    std::make_pair("air_cp", "1006"),
    std::make_pair("hvac_restore_time", "60"),
    std::make_pair("interior_air_min_temp", "293.15"),
    std::make_pair("interior_air_max_temp", "297.039"),
    std::make_pair("air_leakage_ach50", "7"),
    std::make_pair("disable_hvac", "false"),

    std::make_pair("floor_roughness_index", "0"),
    std::make_pair("window_shgc", "0.9"),
    std::make_pair("floor_absorptivity", "0.65"),
    std::make_pair("window_shading_angle", "-1"),  // if the angle between the sun and the z axis is greater than
                                                   // this angle, no direct normal radiation enteres the windows
                                                   // set a negative number for no shading
    std::make_pair("window_areas", "[1.114836, 1.114836, 1.114836, 1.114836]"),
    std::make_pair("window_r_value", "0.5283305514"),

    std::make_pair("house_xlength", "9.23"),
    std::make_pair("house_ylength", "15.38"),
    std::make_pair("house_zlength", "2.46"),
    std::make_pair("house_num_el_x", "5"),
    std::make_pair("house_num_el_y", "8"),
    std::make_pair("house_num_el_z", "4"),

      // 3.5 inches of insulation
    std::make_pair("house_horizontal_thicknesses", "[0.0897]"),
    std::make_pair("house_horizontal_numels",      "[5]"),
    std::make_pair("house_horizontal_kappas",      "[0.039]"),
    std::make_pair("house_horizontal_rhos",        "[45]"),
    std::make_pair("house_horizontal_Cps",         "[2020]"),

    // 6 inches of insulation
    std::make_pair("house_ceiling_thicknesses", "[0.1538]"),
    std::make_pair("house_ceiling_numels",      "[5]"),
    std::make_pair("house_ceiling_kappas",      "[0.039]"),
    std::make_pair("house_ceiling_rhos",        "[45]"),
    std::make_pair("house_ceiling_Cps",         "[2020]"),  

    // 6 inches of concrete
    std::make_pair("house_foundation_thicknesses", "[0.1538]"),
    std::make_pair("house_foundation_numels",      "[5]"),
    std::make_pair("house_foundation_kappas",      "[2.25]"),
    std::make_pair("house_foundation_rhos",        "[2400]"),
    std::make_pair("house_foundation_Cps",         "[880]"),  

    std::make_pair("house_foundation_insulation_thicknesses", "[0.0897]"),
    std::make_pair("house_foundation_insulation_numels",      "[5]"),
    std::make_pair("house_foundation_insulation_kappas",      "[2.1]"),
    std::make_pair("house_foundation_insulation_rhos",        "[2200]"),
    std::make_pair("house_foundation_insulation_Cps",         "[1144]"), 

    std::make_pair("house_ground_horizontal_thickness",  "5"),
    std::make_pair("house_ground_horizontal_numel",      "10"),

    std::make_pair("house_ground_depth",            "20"),
    std::make_pair("house_ground_depth_numel",      "10"),

    // Values from https://open.library.okstate.edu/rainorshine/chapter/13-2-soil-thermal-properties/
    // Table 13-1, for a 2/3 1/3 mixture of clay and soil organic matter  
    std::make_pair("house_ground_kappa",  "2.1"),
    std::make_pair("house_ground_rho",    "2200"),
    std::make_pair("house_ground_Cp",     "1144"),

    std::make_pair("t_start", "0"),
    std::make_pair("t_end", "86400"),
    std::make_pair("timestep_values", "[300, 1800, 1800]"),
    std::make_pair("timestep_points", "[0, 86400, 900000]"),
    std::make_pair("nonlinear_abs_tol", "1e-9"),
    std::make_pair("nonlinear_rel_tol", "1e-8"),
    std::make_pair("nonlinear_itermax", "30"),
    std::make_pair("linear_abs_tol", "1e-12"),
    std::make_pair("linear_rel_tol", "1e-50"),
    std::make_pair("vis_output_freq", "-1"),

    // solar thermal
    std::make_pair("solar_collector_area", "0.0"),
    std::make_pair("solar_collector_efficiency", "0.8"),
    std::make_pair("solar_collector_emissivity", "0.8"),
    std::make_pair("solar_collector_normal", "[0, 0, 1]"),
    std::make_pair("solar_min_thickness", "0.0254"),
  };

  return input_file_defaults;
}

void checkInRange(double val, double val_min, double val_max, const std::string& msg)
{
  if (val < val_min || val > val_max)
    throw std::runtime_error(msg + ", val = " + std::to_string(val) + " must be in range " + std::to_string(val_min) + ", " + std::to_string(val_max));
}

void validateParams(Params& params)
{
  double max = std::numeric_limits<double>::max();

  if (params.weather_filename.size() == 0)
    throw std::runtime_error("weather file must be provided");

  checkInRange(params.ext_wall_emittance, 0, 1, "exterior wall emittance out of range");
  checkInRange(params.ext_wall_absorptivity, 0, 1, "exterior wall absorptivity out of range");

  checkInRange(params.roof_emittance, 0, 1, "roof emittance out of range");
  checkInRange(params.roof_absorptivity, 0, 1, "roof absorptivity out of range");

  checkInRange(params.lawn_emittance, 0, 1, "lawn emittance out of range");
  checkInRange(params.lawn_absorptivity, 0, 1, "lawn absorptivity out of range");

  checkInRange(params.air_rho, 0, max, "air rho out of range");
  checkInRange(params.air_cp, 0, max, "air Cp out of range");
  checkInRange(params.hvac_restore_time, 0, max, "hvac restore time out of range");
  checkInRange(params.interior_air_min_temp, 0, max, "interior air min temp out of range");
  checkInRange(params.interior_air_max_temp, params.interior_air_min_temp, max, "interior air max temp out of range");
  checkInRange(params.interior_air_min_temp, 0, params.interior_air_max_temp, "interior air max temp out of range");
  checkInRange(params.air_leakage_ach50, 0, max, "air leakage ach50 out of range");

  checkInRange(params.window_shgc, 0, max, "window shgc out of range");
  checkInRange(params.floor_absorptivity, 0, 1, "floor absorptivity out of range");
  for (int i=0; i < 4; ++i)
    checkInRange(params.window_areas[i], 0, max, "air rho out of range");

  checkInRange(params.window_r_value, 0, max, "window r value out of range");

}


SimpleHouseSpecSection parseSimpleHouseSection(const std::map<std::string, std::string>& input_vals, const std::string& prefix)
{
  ValueParser parser;
  SimpleHouseSpecSection spec;

  spec.thicknesses = parser.parseArray<double>(input_vals.at(prefix + "_thicknesses"));
  spec.numels      = parser.parseArray<int>(input_vals.at(prefix + "_numels"));
  std::vector<double> kappas      = parser.parseArray<double>(input_vals.at(prefix + "_kappas"));
  std::vector<double> rhos        = parser.parseArray<double>(input_vals.at(prefix + "_rhos"));
  std::vector<double> Cps         = parser.parseArray<double>(input_vals.at(prefix + "_Cps"));

  size_t num_layers = spec.thicknesses.size();
  assertAlways(spec.numels.size() == num_layers, prefix + "_numels length is incorrect");
  assertAlways(kappas.size() == num_layers, "horizontal_kappas length is incorrect");
  assertAlways(rhos.size() == num_layers, "horizontal_rhos length is incorrect");
  assertAlways(Cps.size() == num_layers, "horizontal_Cps length is incorrect");

  for (size_t i=0; i < num_layers; ++i)
    spec.params.emplace_back(kappas[i], rhos[i], Cps[i]);

  return spec;
}


SimpleHouseSpec parseSimpleHouseSpec(const std::map<std::string, std::string>& input_vals)
{
  ValueParser parser;
  SimpleHouseSpec spec;
  double xmax = parser.parseScalar<double>(input_vals.at("house_xlength"));
  double ymax = parser.parseScalar<double>(input_vals.at("house_ylength"));
  double zmax = parser.parseScalar<double>(input_vals.at("house_zlength"));

  int numel_x = parser.parseScalar<int>(input_vals.at("house_num_el_x"));
  int numel_y = parser.parseScalar<int>(input_vals.at("house_num_el_y"));
  int numel_z = parser.parseScalar<int>(input_vals.at("house_num_el_z"));
  spec.middle_block = Mesh::getMeshSpec(0, xmax, 0, ymax, 0, zmax, numel_x, numel_y, numel_z);

  SimpleHouseSpecSection section_spec = parseSimpleHouseSection(input_vals, "house_horizontal");
  spec.horizontal_thicknesses = section_spec.thicknesses;
  spec.horizontal_numels = section_spec.numels;
  spec.horizontal_params = section_spec.params;

  section_spec = parseSimpleHouseSection(input_vals, "house_ceiling");
  spec.ceiling_thicknesses = section_spec.thicknesses;
  spec.ceiling_numels = section_spec.numels;
  spec.ceiling_params = section_spec.params;

  section_spec = parseSimpleHouseSection(input_vals, "house_foundation");
  spec.foundation_thicknesses = section_spec.thicknesses;
  spec.foundation_numels = section_spec.numels;
  spec.foundation_params = section_spec.params;  

  section_spec = parseSimpleHouseSection(input_vals, "house_foundation_insulation");
  assertAlways(section_spec.numels.size() == 1, "Only 1 layer of foundation is allowed");
  spec.foundation_insulation_thicknesses = section_spec.thicknesses;
  spec.foundation_insulation_numels = section_spec.numels;
  spec.foundation_insulation_params = section_spec.params[0];

  spec.ground_horizontal_thickness = parser.parseScalar<double>(input_vals.at("house_ground_horizontal_thickness"));
  spec.ground_horizontal_numel     = parser.parseScalar<int>(input_vals.at("house_ground_horizontal_numel"));

  spec.ground_depth       = parser.parseScalar<double>(input_vals.at("house_ground_depth"));
  spec.ground_depth_numel = parser.parseScalar<int>(input_vals.at("house_ground_depth_numel"));

  double ground_kappa = parser.parseScalar<double>(input_vals.at("house_ground_kappa"));
  double ground_rho = parser.parseScalar<double>(input_vals.at("house_ground_rho"));
  double ground_Cp = parser.parseScalar<double>(input_vals.at("house_ground_Cp"));

  spec.ground_params = Heat::VolumeGroupParams(ground_kappa, ground_rho, ground_Cp);

  return spec;
}

timesolvers::TimeStepperOpts parseTimesolverData(const std::map<std::string, std::string>& input_vals)
{
  ValueParser parser;
  timesolvers::TimeStepperOpts opts;

  opts.t_start = parser.parseScalar<double>(input_vals.at("t_start"));
  opts.t_end   = parser.parseScalar<double>(input_vals.at("t_end"));
  opts.nonlinear_abs_tol = parser.parseScalar<double>(input_vals.at("nonlinear_abs_tol"));
  opts.nonlinear_rel_tol = parser.parseScalar<double>(input_vals.at("nonlinear_rel_tol"));
  opts.nonlinear_itermax = parser.parseScalar<int>(input_vals.at("nonlinear_itermax"));
  opts.vis_output_freq   = parser.parseScalar<int>(input_vals.at("vis_output_freq"));

  std::vector<double> timestep_vals = parser.parseArray<double>(input_vals.at("timestep_values"));
  std::vector<double> timestep_points = parser.parseArray<double>(input_vals.at("timestep_points"));
  assertAlways(timestep_vals.size() == timestep_points.size(), "number of timestep values and points must be the same");
  std::vector<timesolvers::TimestepControllerPiecewise::TimestepPoint> timesteps;
  for (size_t i=0; i < timestep_vals.size(); ++i)
    timesteps.push_back({timestep_points[i], timestep_vals[i]});

  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerPiecewise>(timesteps, commRank(MPI_COMM_WORLD) == 0);
  opts.mat_type = linear_system::LargeMatrixType::Petsc;


  auto matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>();
  matrix_opts->is_structurally_symmetric = true;
  opts.solve_auxiliary_equations_combined_system = true;
  opts.precompute_linear_jacobian = true;

  matrix_opts->petsc_opts["ksp_atol"] = input_vals.at("linear_abs_tol");
  matrix_opts->petsc_opts["ksp_rtol"] = input_vals.at("linear_rel_tol");
  matrix_opts->petsc_opts["ksp_monitor"] = "";

  matrix_opts->petsc_opts["ksp_type"] = "cg";

  if (commSize(MPI_COMM_WORLD) > 1)
  {
    matrix_opts->petsc_opts["ksp_type"] = "cg";
    //matrix_opts->petsc_opts["pc_type"] = "jacobi";

    matrix_opts->petsc_opts["pc_type"] = "asm";
    matrix_opts->petsc_opts["pc_asm_overlap"] = "1";
  }

  opts.matrix_opts = matrix_opts;

  opts.t_start = parser.parseScalar<double>(input_vals.at("t_start"));

  return opts;
}


Params parseParams(const std::string& fname)
{
  InputParser input_parser(fname);
  std::map<std::string, std::string> input_vals = input_parser.parse(getInputFileDefaults());

  Params params;
  ValueParser parser;

  const double pi = std::atan(1)*4;

  params.weather_filename         = input_vals.at("weather_filename");
  params.met_terrain_index        = parser.parseScalar<int>(input_vals.at("met_terrain_index"));
  params.meterological_altitude   = parser.parseScalar<int>(input_vals.at("meterological_altitude"));
  params.local_terrain_index      = parser.parseScalar<int>(input_vals.at("local_terrain_index"));

  params.ext_wall_roughness_index = parser.parseScalar<int>(input_vals.at(   "ext_wall_roughness_index"));
  params.ext_wall_emittance       = parser.parseScalar<double>(input_vals.at("ext_wall_emittance"));
  params.ext_wall_absorptivity    = parser.parseScalar<double>(input_vals.at("ext_wall_absorptivity"));

  params.roof_roughness_index     = parser.parseScalar<int>(input_vals.at(   "roof_roughness_index"));
  params.roof_emittance           = parser.parseScalar<double>(input_vals.at("roof_emittance"));
  params.roof_absorptivity        = parser.parseScalar<double>(input_vals.at("roof_absorptivity"));

  params.lawn_roughness_index     = parser.parseScalar<int>(input_vals.at(   "lawn_roughness_index"));
  params.lawn_emittance           = parser.parseScalar<double>(input_vals.at("lawn_emittance"));
  params.lawn_absorptivity        = parser.parseScalar<double>(input_vals.at("lawn_absorptivity"));

  params.int_wall_roughness_index = parser.parseScalar<int>(input_vals.at("int_wall_roughness_index"));

  params.air_rho                  = parser.parseScalar<double>(input_vals.at("air_rho"));
  params.air_cp                   = parser.parseScalar<double>(input_vals.at("air_cp"));
  params.hvac_restore_time        = parser.parseScalar<double>(input_vals.at("hvac_restore_time"));
  params.interior_air_min_temp    = parser.parseScalar<double>(input_vals.at("interior_air_min_temp"));
  params.interior_air_max_temp    = parser.parseScalar<double>(input_vals.at("interior_air_max_temp"));
  params.air_leakage_ach50        = parser.parseScalar<double>(input_vals.at("air_leakage_ach50"));
  params.disable_hvac             = parser.parseScalar<bool>(input_vals.at("disable_hvac"));


  params.floor_roughness_index     = parser.parseScalar<int>(input_vals.at(   "floor_roughness_index"));
  params.window_shgc               = parser.parseScalar<double>(input_vals.at("window_shgc"));
  params.floor_absorptivity        = parser.parseScalar<double>(input_vals.at("floor_absorptivity"));
  params.window_shading_angle      = parser.parseScalar<double>(input_vals.at("window_shading_angle"));
  params.window_shading_angle      = params.window_shading_angle * pi / 180.0;
  std::vector<double> window_areas = parser.parseArray<double>(input_vals.at( "window_areas"));
  if (window_areas.size() != 4)
    throw std::runtime_error("the building has 4 walls, must have exactly 4 window areas");

  params.window_areas              = {window_areas[0], window_areas[1], window_areas[2], window_areas[3]};
  params.window_area               = window_areas[0] + window_areas[1] + window_areas[2] + window_areas[3];
  params.window_r_value            = parser.parseScalar<double>(input_vals.at("window_r_value"));

  
  params.solar_collector_area       = parser.parseScalar<double>(input_vals.at(  "solar_collector_area"));
  params.solar_collector_efficiency = parser.parseScalar<double>(input_vals.at(  "solar_collector_efficiency"));
  params.solar_collector_emissivity = parser.parseScalar<double>(input_vals.at(  "solar_collector_emissivity"));
  auto solar_collector_normal       = parser.parseArray<double>(input_vals.at("solar_collector_normal"));
  assertAlways(solar_collector_normal.size() == 3, "solar_collect_normal must have 3 components");
  params.solar_collector_normal     = {solar_collector_normal[0], solar_collector_normal[1], solar_collector_normal[2]};
  params.solar_collector_normal     = params.solar_collector_normal / std::sqrt(dot(params.solar_collector_normal, params.solar_collector_normal));
  params.solar_min_thickness        = parser.parseScalar<double>(input_vals.at(  "solar_min_thickness"));


  params.simple_house_spec = parseSimpleHouseSpec(input_vals);
  params.time_stepper_opts = parseTimesolverData(input_vals);

  validateParams(params);

  return params;
}

}  // namespace