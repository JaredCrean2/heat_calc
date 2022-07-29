#ifndef PHYSICS_HEAT_INTERIOR_LOADS_H
#define PHYSICS_HEAT_INTERIOR_LOADS_H

#include "ProjectDefs.h"
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



}  // namespace

#endif