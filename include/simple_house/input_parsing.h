#ifndef HEAT_CALC_SIMPLE_HOUSE_INPUT_PARSING_H
#define HEAT_CALC_SIMPLE_HOUSE_INPUT_PARSING_H

#include "ProjectDefs.h"
#include "file/input_parser.h"
#include "linear_system/large_matrix_petsc.h"
#include "simple_house_spec.h"
#include "time_solver/time_stepper_opts.h"
#include "time_solver/timestep_controller_piecewise.h"

namespace simple_house {


struct Params
{
  std::string weather_filename;

  // for exterior wall  + lawn Tarp BCs
  std::array<Real, 3> vertical_vector = {0, 0, 1};
  std::array<Real, 3> point_at_zero_altitude = {0, 0, -100};  
  int met_terrain_index = 2;  // rough, wooded country
  Real meterological_altitude = 1836;
  int local_terrain_index = 2;  // rough, wooded country

  // for exterior wall BCs
  // absorptivity and emissivity values for stucco from: https://remdb.nrel.gov/measures.php?gId=12&ctId=216&scId=2374
  // absorptivity and emissivity values for asphalt singles from Medina "Effects of Single Absorptivity, 
  // radient barrier emissivity", International Journal of Energy Research, 2000, 24:665  
  int ext_wall_roughness_index = 0;  //TODO: check this value
  double ext_wall_emittance = 0.9;
  double ext_wall_absorptivity = 0.75;

  // for roof BC
  int roof_roughness_index = 0; // TODO: check this value
  Real roof_emittance = 0.78;
  Real roof_absorptivity = 0.78;

  // for lawn BC
  int lawn_roughness_index = 0;
  Real lawn_emittance    = 0.95;  // Van Wijk and Scholte Ubing (1963)
  Real lawn_absorptivity = 0.8;   // could not find a great source for this, different forms of clay
                                  // have absorptivities in the 0.7 - 0.9 range  

  // for interior wall
  int int_wall_roughness_index = 0;

  // for HVAC
  Real air_rho           = 1.007;
  Real air_cp            = 1006;
  Real hvac_restore_time = 60;  
  Real interior_air_min_temp = 293.15;
  Real interior_air_max_temp = 297.039; 
  Real air_leakage_ach50 = 7;

  // for windows + floor radiation BC
  int floor_roughness_index = 0;
  Real window_shgc = 0.9;
  Real floor_absorptivity = 0.65; // unfinished concrete 
  Real window_shading_angle = -1;
  std::array<Real, 4> window_areas = {0.557418 * 2, 0.557418 * 2, 0.557418 * 2, 0.557418 * 2};  // 6 sq ft each window, 2 windows per wall
  Real window_area    = window_areas[0] + window_areas[1] + window_areas[2] + window_areas[3];
  Real window_r_value = 3 * 0.1761101838;  // r value converted to SI units   

  SimpleHouseSpec simple_house_spec;
  timesolvers::TimeStepperOpts time_stepper_opts;                        
};

const std::map<std::string, std::string>& getInputFileDefaults();

void checkInRange(double val, double val_min, double val_max, const std::string& msg);

void validateParams(const Params& params);

struct SimpleHouseSpecSection
{
  std::vector<double> thicknesses;
  std::vector<int> numels;
  std::vector<Heat::VolumeGroupParams> params;
};

SimpleHouseSpecSection parseSimpleHouseSection(const std::map<std::string, std::string>& input_vals, const std::string& prefix);

SimpleHouseSpec parseSimpleHouseSpec(const std::map<std::string, std::string>& input_vals);

timesolvers::TimeStepperOpts parseTimesolverData(const std::map<std::string, std::string>& input_vals);

Params parseParams(const std::string& fname);

}

#endif