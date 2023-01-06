#ifndef PHYSICS_HEAT_SOLAR_POSITION
#define PHYSICS_HEAT_SOLAR_POSITION

#include "ProjectDefs.h"
#include "utils/math.h"
#include "dates.h"

namespace Heat {

// azimuth is angle in the E-N coordinate system, measured from north
// Note: East is the x direction and North is the Y direction
// zenith angle from vertical direction to the sun
struct AzimuthZenith
{
  Real cos_azimuth;
  Real cos_zenith;
};


struct DirectionCosines
{
  Real cs1;
  Real cs2;
  Real cs3;
};

namespace solar {

// longitude: radians, east is positive
// time_zone: hours from UTC (ex. Mountain time is UTC - 7, time_zone = 7)
AzimuthZenith computeAzimuthZenith_lowprecision(const Date& date, Real hour, int time_zone, Real latitude, Real longitude);


struct DecTimeDist
{
  Real declination_angle = 0;  // delta in Walton
  Real equation_of_time  = 0;  // epsilon in Walton
  Real distance_to_sun   = 0;  // R in Walton
};

DecTimeDist computeDecTimeDist(Real julian_day);

// time_hours is the time in hours since midnight
// time_zone is the integer identifying the time zone.
//   Ex. Eastern time is GMT - 5.  time_zone is 5
// longtiude: longtiude (in radians), west is positive
Real computeHourAngle(Real time_hours, const DecTimeDist& dec_time_dist, int time_zone, Real longitude);

// latitude: west is positive
DirectionCosines computeDirectionCosines(const DecTimeDist& dec_time_dist, Real hour_angle, Real latitude);

DirectionCosines computeDirectionCosines(int julian_day, int time_hours, int time_zone, Real latitude, Real longitude);

DirectionCosines computeDirectionCosines(const Date& date, int time_hours, int time_zone, Real latitude, Real longitude);

AzimuthZenith computeAzimuthZenith(const DirectionCosines& cosines);

AzimuthZenith computeAzimuthZenith(int julian_day, int time_hours, int time_zone, Real latitude, Real longitude);

AzimuthZenith computeAzimuthZenith(const Date& date, int time_hours, int time_zone, Real latitude, Real longitude);


// converts angles specified in Degrees, Minutes, Seconds to decimal (radians)
Real DMSToRadians(Real degrees, Real minutes, Real seconds);

}  // namespace solar

}  // namespace Heat

#endif