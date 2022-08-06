#ifndef PHYSICS_HEAT_SOLAR_POSITION_CALCULATOR_H
#define PHYSICS_HEAT_SOLAR_POSITION_CALCULATOR_H

#include "ProjectDefs.h"
#include "solar_position.h"

namespace Heat {

class SolarPositionCalculator
{
  public:
    SolarPositionCalculator(Real julian_date_start, int time_zone, Real longitude, Real latitude) :
    m_julian_date_start(julian_date_start),
      m_time_zone(time_zone),
      m_longitude(longitude),
      m_latitude(latitude)
    {}

    DirectionCosines computePosition(int t)
    {
      Real julian_date = m_julian_date_start + t;
      int julian_date_whole = static_cast<int>(julian_date);
      Real time_hours = julian_date_whole - julian_date;
      return solar::computeDirectionCosines(julian_date_whole, time_hours, m_time_zone, m_longitude, m_latitude);
    }

  private:
    Real m_julian_date_start;
    int m_time_zone;
    Real m_longitude;
    Real m_latitude;
};

}

#endif