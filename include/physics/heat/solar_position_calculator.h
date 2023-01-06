#ifndef PHYSICS_HEAT_SOLAR_POSITION_CALCULATOR_H
#define PHYSICS_HEAT_SOLAR_POSITION_CALCULATOR_H

#include "ProjectDefs.h"
#include "solar_position.h"

#include <iostream>

namespace Heat {

class SolarPositionCalculator
{
  public:
    // time_zone is the integer identifying the time zone.
    //   Ex. Eastern time is GMT - 5.  time_zone is 5
   // longtiude and latitude are in radians, north and west are positive 
    SolarPositionCalculator(Real julian_date_start, int time_zone, Real latitude, Real longitude) :
    m_julian_date_start(julian_date_start),
      m_time_zone(time_zone),
      m_latitude(latitude),
      m_longitude(longitude)
    {}

    DirectionCosines computePositionFromSeconds(Real t_seconds)
    {
      auto cosines = computePosition(t_seconds/(60*60*24));
      std::cout << "solar angles = " << std::acos(cosines.cs1)*180/3.14159265 << ", " << std::acos(cosines.cs2)*180/3.14159265 << ", " << std::acos(cosines.cs3)*180/3.14159265 << std::endl;
      std::cout << "direction cosines = " << cosines.cs1 << ", " << cosines.cs1 << ", " << cosines.cs3 << std::endl;
      return cosines;
    }

    // t is in days
    DirectionCosines computePosition(Real t)
    {
      std::cout << "computing solar position at t = " << t << std::endl;
      Real julian_date = m_julian_date_start + static_cast<Real>(t);
      std::cout << "julian_date = " << julian_date << std::endl;
      int julian_date_whole = static_cast<int>(julian_date);
      Real time_hours = 24*(julian_date_whole - julian_date);
      return solar::computeDirectionCosines(julian_date_whole, time_hours, m_time_zone, m_latitude, m_longitude);
    }

  private:
    Real m_julian_date_start;
    int m_time_zone;
    Real m_latitude;
    Real m_longitude;
};

}

#endif