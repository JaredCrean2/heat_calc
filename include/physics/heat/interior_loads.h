#ifndef PHYSICS_HEAT_INTERIOR_LOADS_H
#define PHYSICS_HEAT_INTERIOR_LOADS_H

#include "ProjectDefs.h"
#include "physics/post_processor_base.h"

// heat sources within the building (that increase air temperature)

namespace Heat {

class InteriorLoads
{
  public:
    virtual ~InteriorLoads() {};

    // Power of all interior loads (W).  Positive means interior air temperature
    // increases.
    virtual Real computeLoadPower() = 0;
};

class InteriorLoadsConstant : public InteriorLoads
{
  public:
    InteriorLoadsConstant(Real load) :
      m_load(load)
    {}

    Real computeLoadPower() override { return m_load; }
  
  private:
    Real m_load;
};

class PostProcessorInteriorLoads : public physics::PostProcessorBase
{
  public:
    PostProcessorInteriorLoads(std::shared_ptr<InteriorLoads> model) :
      m_model(model)
    {}

    // returns number of values this postprocessor returns
    virtual int numValues() const { return 1; }

    virtual std::vector<std::string> getNames() const { return {"interior_loads"}; }

    virtual std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
    {
      return {m_model->computeLoadPower()};
    }

  private:
    std::shared_ptr<InteriorLoads> m_model;
};



}  // namespace

#endif