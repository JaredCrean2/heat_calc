#ifndef PHYSICS_HEAT_SOLAR_POSITION_CALCULATOR_H
#define PHYSICS_HEAT_SOLAR_POSITION_CALCULATOR_H

#include "ProjectDefs.h"
#include "solar_position.h"

#include <iostream>

namespace Heat {

class SolarPositionCalculator
{
  public:
    virtual ~SolarPositionCalculator() {}

    virtual DirectionCosines computePosition(Real t_seconds) = 0;

};


class SolarPositionCalculatorConstant : public SolarPositionCalculator
{
  public:
    SolarPositionCalculatorConstant(const DirectionCosines& cosines) :
      m_cosines(cosines)
    {}

    DirectionCosines computePosition(Real t_seconds) override
    {
      return m_cosines;
    }



  private:
    DirectionCosines m_cosines;
};

// based on the Naval Observatory algorithm
class SolarPositionCalculatorNaval : public SolarPositionCalculator
{
  public:
    // time_zone is the integer identifying the time zone.
    //   Ex. Eastern time is GMT - 5.  time_zone is 5
   // longtiude and latitude are in radians, north and west are positive 
    SolarPositionCalculatorNaval(Real julian_date_start, int time_zone, Real latitude, Real longitude) :
    m_julian_date_start(julian_date_start),
      m_time_zone(time_zone),
      m_latitude(latitude),
      m_longitude(longitude)
    {}

    // t is in seconds from the starting time
    DirectionCosines computePosition(Real t_seconds) override
    {
      auto cosines = computePositionFromDays(t_seconds/(60*60*24));

      return cosines;
    }

  private:

    // t is in days
    DirectionCosines computePositionFromDays(Real t)
    {
      Real julian_date = m_julian_date_start + static_cast<Real>(t);
      int julian_date_whole = static_cast<int>(julian_date);
      Real time_hours = 24*(julian_date_whole - julian_date);
      return solar::computeDirectionCosines(julian_date_whole, time_hours, m_time_zone, m_latitude, m_longitude);
    }

    Real m_julian_date_start;
    int m_time_zone;
    Real m_latitude;
    Real m_longitude;
};

}

#endif