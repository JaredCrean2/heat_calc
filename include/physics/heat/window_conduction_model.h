#ifndef PHYSICS_HEAT_WINDOW_CONDUCTION_MODEL_H
#define PHYSICS_HEAT_WINDOW_CONDUCTION_MODEL_H

#include "ProjectDefs.h"

namespace Heat {

class WindowConductionModel
{
  public:
    WindowConductionModel(Real r_val, Real area) :
      m_r_val(r_val),
      m_area(area)
    {}

    // computes rate of energy loss (W) via conduction through windows.
    // Positive values mean energy is adding to the interior air
    Real computeConductionPower(Real t_interior, Real t_exterior)
    {
      return (t_exterior - t_interior) * m_area/m_r_val;
    }


  private:
    Real m_r_val;  // (m^2 * K)/W
    Real m_area;
};

}  // namespace


#endif