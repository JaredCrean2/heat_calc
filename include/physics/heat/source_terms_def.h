#ifndef PHYSICS_HEAT_SOURCE_TERM_DEFS_H
#define PHYSICS_HEAT_SOURCE_TERM_DEFS_H

#include "discretization/source_term.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_radiation.h"
#include "utils/bounding_box.h"


namespace Heat {


class SourceTermAirWindSky : public SourceTerm
{
  public:
    SourceTermAirWindSky(VolDiscPtr vol_disc) :
      SourceTerm(vol_disc)
    {}

    virtual void setAirTemperature(Real temp) {}
    
    virtual void setAirSpeed(Real velocity) {}
    
    virtual void setAirDirection(std::array<Real, 3> direction) {}

    virtual void setIRHorizontalRadiation(Real flux) {}

    virtual void setDirectNormalRadiation(Real flux) {}

    virtual void setDiffuseRadiation(Real flux) {}

    virtual void setSolarDirection(const Heat::DirectionCosines& cosines) {}  
};

constexpr Real MIN_REAL = std::numeric_limits<Real>::lowest();
constexpr Real MAX_REAL = std::numeric_limits<Real>::max();


// Models solar thermal heating by modeling a solar collector and 
// distributing the flux evenly over a given area
class SourceTermSolarHeating : public SourceTermAirWindSky
{
  public:

    SourceTermSolarHeating(VolDiscPtr vol_disc, Real collector_area, Real collector_efficiency,
                           const std::array<Real, 3>& collector_unit_normal,
                           Real output_area, const utils::BoundingBox& box) :
      SourceTermAirWindSky(vol_disc),
      m_solar_model(collector_efficiency),
      m_collector_area(collector_area),
      m_unit_normal(collector_unit_normal),
      m_output_area(output_area),
      m_box(box)
    {}
    
    SourceTermSolarHeating(VolDiscPtr vol_disc, Real collector_area, Real collector_efficiency,
                           const std::array<Real, 3>& collector_unit_normal,
                           Real output_area) :
      SourceTermSolarHeating(vol_disc, collector_area, collector_efficiency, collector_unit_normal,
                             output_area, utils::BoundingBox({MIN_REAL, MIN_REAL, MIN_REAL}, {MAX_REAL, MAX_REAL, MAX_REAL}))
    {}
    
    void setDirectNormalRadiation(Real flux) override { m_solar_model.setDirectNormalRadiation(flux); }

    void setDiffuseRadiation(Real flux) override { m_solar_model.setDiffuseRadiation(flux); }

    void setSolarDirection(const Heat::DirectionCosines& cosines) override { m_solar_model.setSolarDirection(cosines); }  

    void getValues(Index el, const Real t, Real* vals) override
    {
      Real flux = m_solar_model.computeFlux(m_unit_normal);
      flux = flux * m_collector_area / m_output_area;
      bool in_box = m_box.contains(computeCentroid(el));

      for (int i=0; i < getVolDisc()->getNumQuadPtsPerElement(); ++i)
        vals[i] = in_box ? flux : 0;
    }


  private:

    std::array<Real, 3> computeCentroid(Index el)
    {
      std::array<Real, 3> centroid = {0, 0, 0};
      for (int i=0; i < getVolDisc()->getNumCoordPtsPerElement(); ++i)
        for (int d=0; d < 3; ++d)
          centroid[d] += getVolDisc()->vol_group.coords[el][i][d];

      return centroid;
    }
    Heat::SolarRadiationModel m_solar_model;
    Real m_collector_area;
    std::array<Real, 3> m_unit_normal;
    Real m_output_area;
    utils::BoundingBox m_box;
};

}

#endif