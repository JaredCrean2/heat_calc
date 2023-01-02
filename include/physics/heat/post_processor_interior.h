#ifndef PHYSISC_HEAT_POST_PROCESSOR_INTERIOR_H
#define PHYSISC_HEAT_POST_PROCESSOR_INTERIOR_H

#include "physics/heat/HeatEquationSolar.h"
#include "physics/post_processor_base.h"

namespace Heat {

class InteriorAirTemperatureUpdator;

class AuxiliaryEquationsSolar;
using AuxiliaryEquationsSolarPtr = std::shared_ptr<AuxiliaryEquationsSolar>;

class PostProcessorInterior : public physics::PostProcessorBase
{
  public:
    PostProcessorInterior(AuxiliaryEquationsSolarPtr aux_eqns, std::shared_ptr<InteriorAirTemperatureUpdator> air_temp) :
      m_aux_eqns(aux_eqns),
      m_air_temp(air_temp)
    {}

    int numValues() const override { return 2; }

    std::vector<std::string> getNames() const override { return {"interior_air_temp", "hvac_flux"}; }

    std::vector<Real> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override;

  private:
    AuxiliaryEquationsSolarPtr m_aux_eqns;
    std::shared_ptr<InteriorAirTemperatureUpdator> m_air_temp;
};

}  // namespace

#endif