#include "physics/heat/solar_position.h"
#include "error_handling.h"
#include "math.h"
#include "physics/heat/dates.h"
#include <iostream>

namespace Heat {

namespace solar {


AzimuthZenith computeAzimuthZenith_lowprecision(const Date& date, Real hour, int time_zone, Real latitude, Real longitude)
{
  // based on https://gml.noaa.gov/grad/solcalc/solareqns.PDF
  assertAlways(longitude >= -PI && longitude <= PI, "longitude must be in range [-pi, pi] radians");

  // conver to the conventions used by the algorithms
  time_zone = -time_zone;
  longitude = -longitude;


  int days_per_year = isLeapYear(date.year) ? 366 : 365;
  int day_of_year = computeDayOfYear(date);

  Real day_fraction = (hour - 12)/24;
  Real gamma = 2*PI*(day_of_year - 1 + day_fraction)/days_per_year;

  Real equation_of_time = 229.18*(0.000075 + 0.001868*std::cos(gamma) - 0.032077*std::sin(gamma) - 
                          0.014615*std::cos(2*gamma) - 0.040849*std::sin(2*gamma));

  Real declination_angle = 0.006918 - 0.399912*std::cos(gamma) + 0.070257*std::sin(gamma) - 
                           0.006758*std::cos(2*gamma) + 0.000907*std::sin(2*gamma) - 
                           0.002697*std::cos(3*gamma) + 0.00148*std::sin(3*gamma);

  Real time_offset_minutes = equation_of_time + 4*radiansToDegrees(longitude) - 60*time_zone;
  Real tst = hour*60 + time_offset_minutes;  // true solar time?

  Real hour_angle = degreesToRadians((tst/4) - 180);

  Real cos_zenith = std::sin(latitude)*std::sin(declination_angle) + std::cos(latitude)*std::cos(declination_angle)*std::cos(hour_angle);
  Real sin_zenith = std::sin(std::acos(cos_zenith));
  Real azimuth = -std::acos(-(std::sin(latitude)*cos_zenith - std::sin(declination_angle))/(std::cos(latitude)*sin_zenith));

  AzimuthZenith az;
  az.cos_azimuth = std::cos(azimuth);
  az.cos_zenith  = cos_zenith;

  return az;
}


DecTimeDist computeDecTimeDist(int julian_date)
{
  Real T = (julian_date - 2415020.5) / 36525;
  Real G = 385.475833 + 35999.049750*T - 0.000150*T*T;   G = degreesToRadians(G);
  Real L = 279.696678 + 36000.768920*T + 0.000303*T*T;   L = degreesToRadians(L);
  Real C = 270.434164 + 481267.883144*T - 0.001133*T*T;  C = degreesToRadians(C);
  Real N = 259.183275 - 1934.142008*T + 0.002078*T*T;    N = degreesToRadians(N);
  Real V = 212.603219 + 58517.803875*T + 0.0001286*T*T;  V = degreesToRadians(V);
  Real M = 319.529425 + 19139.858500*T + 0.000181*T*T;   M = degreesToRadians(M);
  Real J = 225.444651 + 3034.906654*T;                   J = degreesToRadians(J);


  Real theta = 0.397930 * std::sin(L) + 0.009999*std::sin(G-L) + 0.003334 * std::sin(G + L) -
               0.000208*T*std::sin(L) + 0.000042*std::sin(2*G+L) - 0.000040*std::cos(L) -
               0.000039*std::sin(N-L) - 0.000030*T*std::sin(G-L) - 0.000014*std::sin(2*G-L) -
               0.000010*std::cos(G-L-J) - 0.000010*T*std::sin(G+L);

  Real rho = 1.000421 - 0.033503*std::cos(G) - 0.000140*std::cos(2*G) +
             0.000084*T*std::cos(G) - 0.000033*std::sin(G-J) + 0.000027*std::sin(2*G-2*V);

  Real phi = -0.041295*std::sin(2*L) + 0.032116*std::sin(G) - 0.001038*std::sin(G-2*L) -
              0.000346*std::sin(G+2*L) - 0.000095 - 0.000080*T*std::sin(G) -
              0.000079*std::sin(N) + 0.000068*std::sin(2*G) + 0.000046*T*std::sin(2*L)
              + 0.000030*std::sin(C-L) - 0.000025*std::cos(G-J) + 0.000024*std::sin(4*G-8*M+3*J)
              - 0.000019*std::sin(G-V) - 0.000017*std::cos(2*G - 2*V);

  DecTimeDist dec_time_dist;
  dec_time_dist.declination_angle = std::asin(theta / std::sqrt(rho));
  dec_time_dist.equation_of_time  = radiansToDegrees(-std::asin(phi / std::sqrt(rho - theta*theta)))/15;  // unit: hours
  dec_time_dist.distance_to_sun   = std::sqrt(rho);

  return dec_time_dist;
}


Real computeHourAngle(Real time_hours, const DecTimeDist& dec_time_dist, int time_zone, Real longitude)
{
  assertAlways(longitude >= -PI && longitude <= PI, "longitude must be in range [-pi, pi] radians");
  Real longitude_degrees = radiansToDegrees(longitude);
  
  auto val =  15 * (12 - time_hours - dec_time_dist.equation_of_time  - time_zone) + longitude_degrees;
  auto val_radians = degreesToRadians(val);

  return val_radians;
}


DirectionCosines computeDirectionCosines(const DecTimeDist& dec_time_dist, Real hour_angle, Real latitude)
{
  DirectionCosines cosines;
  Real sin_delta = std::sin(dec_time_dist.declination_angle);
  Real cos_delta = std::cos(dec_time_dist.declination_angle);

  Real sin_lambda = std::sin(latitude);
  Real cos_lambda = std::cos(latitude);

  Real sin_h = std::sin(hour_angle);
  Real cos_h = std::cos(hour_angle);

  //Note: the negative sign in front of cosins.cs1 does not appear in the TARP
  //      manual.  It appear the coordinate system in TARP has West as the +x
  //      direction, rather than east.  The negative sign corrects the problem.
  cosines.cs1 = -cos_delta * sin_h;
  cosines.cs2 = sin_delta * cos_lambda - cos_delta * sin_lambda * cos_h;
  cosines.cs3 = sin_delta * sin_lambda + cos_delta * cos_lambda * cos_h;

  return cosines;
}


AzimuthZenith computeAzimuthZenith(const DirectionCosines& cosines)
{
  AzimuthZenith az;
  az.cos_zenith  = cosines.cs3;
  az.cos_azimuth = cosines.cs2 / std::sqrt(1 - cosines.cs3*cosines.cs3);

  return az;
}


DirectionCosines computeDirectionCosines(int julian_day, Real time_hours, int time_zone, Real latitude, Real longitude)
{
  DecTimeDist dec_time_dist = computeDecTimeDist(julian_day);
  Real hour_angle = computeHourAngle(time_hours, dec_time_dist, time_zone, longitude);
  return computeDirectionCosines(dec_time_dist, hour_angle, latitude);
}

DirectionCosines computeDirectionCosines(const Date& date, Real time_hours, int time_zone, Real latitude, Real longitude)
{
  return computeDirectionCosines(computeJulianDate(date), time_hours, time_zone, latitude, longitude);
}

AzimuthZenith computeAzimuthZenith(int julian_day, Real time_hours, int time_zone, Real latitude, Real longitude)
{
  return computeAzimuthZenith(computeDirectionCosines(julian_day, time_hours, time_zone,  latitude, longitude));
}

AzimuthZenith computeAzimuthZenith(const Date& date, Real time_hours, int time_zone, Real latitude, Real longitude)
{
  return computeAzimuthZenith(computeJulianDate(date), time_hours, time_zone, latitude, longitude);
}



Real DMSToRadians(Real degrees, Real minutes, Real seconds)
{
  Real angle_degrees = degrees + minutes/60 +  seconds/3600;
  return degreesToRadians(angle_degrees);
}

} // namespace

}  // namespace