#include "physics/heat/post_processor_environment_data.h"

namespace Heat
{

// returns number of values this postprocessor returns
int PostProcessorEnvironmentData::numValues() const
{ 
  return 9;
}

std::vector<std::string> PostProcessorEnvironmentData::getNames() const 
{ 
  return {"exterior_air_temp", "exterior_air_speed", "air_direction_x", "air_direction_y", "air_direction_z",
          "ir_horizontal_radiation", "direct_normal_radiation", "diffuse_radiation", "total_radiation"}; 
}

std::vector<Real> PostProcessorEnvironmentData::getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
{
  std::vector<Real> vals(numValues());
  auto& edata = m_heat_eqn->getEnvironmentData();

  vals[0] = edata.air_temp;
  vals[1] = edata.air_speed;
  vals[2] = edata.air_direction[0];
  vals[3] = edata.air_direction[1];
  vals[4] = edata.air_direction[2];
  vals[5] = edata.ir_horizontal_radiation;
  vals[6] = edata.direct_normal_radiation;
  vals[7] = edata.diffuse_radiation;
  vals[8] = edata.ir_horizontal_radiation + edata.direct_normal_radiation + edata.diffuse_radiation;

  return vals;
}

}  // namespace