#ifndef PHYSICS_HEAT_POST_PROCESSOR_ENVIRONMENT_DATA_H
#define PHYSICS_HEAT_POST_PROCESSOR_ENVIRONMENT_DATA_H

#include "physics/heat/HeatEquationSolar.h"
#include "physics/post_processor_base.h"

namespace Heat {

class PostProcessorEnvironmentData : public physics::PostProcessorBase
{
  public:
    PostProcessorEnvironmentData(HeatEquationSolar* heat_eqn_solar) :
      m_heat_eqn(heat_eqn_solar)
    {}

    // returns number of values this postprocessor returns
    int numValues() const override;

    std::vector<std::string> getNames() const override;

    std::vector<Real> getValues(DiscVectorPtr u, double t) override;

  private:
    HeatEquationSolar* m_heat_eqn;
};

}  // namespace

#endif