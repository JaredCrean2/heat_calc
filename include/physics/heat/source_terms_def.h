#ifndef PHYSICS_HEAT_SOURCE_TERM_DEFS_H
#define PHYSICS_HEAT_SOURCE_TERM_DEFS_H

#include "discretization/source_term.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_radiation.h"
#include "utils/bounding_box.h"
#include "physics/post_processor_base.h"
#include "physics/heat/temperature_controller.h"


namespace Heat {


class SourceTermAirWindSky : public SourceTerm
{
  public:
    SourceTermAirWindSky(VolDiscPtr vol_disc) :
      SourceTerm(vol_disc)
    {}

    // gets source term derivative wrt air temperature at quadrature nodes
    virtual void getValues_dTair(Index el, const Real t, Real t_air_dot, Real* vals_dot) = 0;

    virtual void setInteriorAirTemperature(Real temp) {}

    virtual void setExteriorAirTemperature(Real temp) {}
    
    virtual void setAirSpeed(Real velocity) {}
    
    virtual void setAirDirection(std::array<Real, 3> direction) {}

    virtual void setIRHorizontalRadiation(Real flux) {}

    virtual void setDirectNormalRadiation(Real flux) {}

    virtual void setDiffuseRadiation(Real flux) {}

    virtual void setSolarDirection(const Heat::DirectionCosines& cosines) {}  
};

using SourceTermAirWindSkyPtr = std::shared_ptr<SourceTermAirWindSky>;


// Models solar thermal heating by modeling a solar collector and 
// distributing the flux evenly over a given area
class SourceTermSolarHeating : public SourceTermAirWindSky
{
  public:

    SourceTermSolarHeating(VolDiscPtr vol_disc, Real collector_area, Real collector_efficiency,
                           Real collector_emissivity,
                           const std::array<Real, 3>& collector_unit_normal,
                           Real output_area, const utils::BoundingBox& box,
                           std::shared_ptr<TemperatureController> controller) :
      SourceTermAirWindSky(vol_disc),
      m_solar_model(1),
      m_collector_efficiency(collector_efficiency),
      m_collector_emissivity(collector_emissivity),
      m_collector_area(collector_area),
      m_unit_normal(collector_unit_normal),
      m_output_area(output_area),
      m_box(box),
      m_controller(controller)
    {  
    }
    
    SourceTermSolarHeating(VolDiscPtr vol_disc, Real collector_area, Real collector_efficiency,
                           Real collector_emissivity,
                           const std::array<Real, 3>& collector_unit_normal,
                           Real output_area) :
      SourceTermSolarHeating(vol_disc, collector_area, collector_efficiency, collector_emissivity, collector_unit_normal,
                             output_area, utils::BoundingBox(),
                             std::make_shared<TemperatureControllerConstant>())
    {
    }
    
    void setIRHorizontalRadiation(Real flux) override { m_t_sky4 = flux/m_sigma; }

    void setInteriorAirTemperature(Real t_air) override
    {
      m_t_air_interior = t_air;
    }

    void setExteriorAirTemperature(Real t_air) override
    {
      m_t_air_exterior = t_air;
      m_t_air_exterior4 = std::pow(t_air, 4);
    }    

    void setDirectNormalRadiation(Real flux) override { m_solar_model.setDirectNormalRadiation(flux); }

    void setDiffuseRadiation(Real flux) override { m_solar_model.setDiffuseRadiation(flux); }

    void setSolarDirection(const Heat::DirectionCosines& cosines) override { m_solar_model.setSolarDirection(cosines); }  

    Real getTotalFlux()
    {
      return m_controller->get_value(m_t_air_interior) * computeCollectorFlux();
    }

    Real getTotalFlux_dTair(Real t_air_dot)
    {
      return m_controller->get_value_dot(m_t_air_interior, t_air_dot) * computeCollectorFlux();
    }  

    Real getControllerFactor()
    {
      return m_controller->get_value(m_t_air_interior);
    }  

    void getValues(Index el, const Real t, Real* vals) override
    {
      Real flux_per_volume = getTotalFlux() / m_output_area;
      bool in_box = m_box.contains(computeCentroid(el));

      for (int i=0; i < getVolDisc()->getNumQuadPtsPerElement(); ++i)
        vals[i] = in_box ? flux_per_volume : 0;
    }

    void getValues_dTair(Index el, const Real t, Real t_air_dot, Real* vals_dot) override
    {
      //std::cout << "t_air_dot = " << t_air_dot << std::endl;
      Real flux_per_volume_dot = getTotalFlux_dTair(t_air_dot) / m_output_area;
      bool in_box = m_box.contains(computeCentroid(el));
      //std::cout << "flux dot = " << flux_per_volume_dot << std::endl;
      //std::cout << std::boolalpha << "in_box = " << in_box << std::endl;

      for (int i=0; i < getVolDisc()->getNumQuadPtsPerElement(); ++i)
        vals_dot[i] = in_box ? flux_per_volume_dot : 0;
    }


  private:
    Real computeCollectorFlux()
    {
      Real flux = m_solar_model.computeFlux(m_unit_normal);  

      // assume that the collector temperature is equal to the air temperature (even
      // though that is only approximate)
      Real sky_flux = m_sigma * m_collector_emissivity * (m_t_sky4 - m_t_air_exterior4);
      Real net_flux = m_collector_efficiency * (flux + sky_flux) * m_collector_area;
      //std::cout << "direct flux = " << flux << ", sky_flux = " << sky_flux << std::endl;

      return std::max(net_flux, 0.0);      
    }

    std::array<Real, 3> computeCentroid(Index el)
    {
      std::array<Real, 3> centroid = {0, 0, 0};
      for (int i=0; i < getVolDisc()->getNumCoordPtsPerElement(); ++i)
        for (int d=0; d < 3; ++d)
          centroid[d] += getVolDisc()->vol_group.coords[el][i][d];

      for (int d=0; d < 3; ++d)
        centroid[d] /= getVolDisc()->getNumCoordPtsPerElement();

      return centroid;
    }

    Heat::SolarRadiationModel m_solar_model;
    Real m_collector_efficiency;
    Real m_collector_emissivity;
    Real m_collector_area;
    std::array<Real, 3> m_unit_normal;
    Real m_output_area;
    utils::BoundingBox m_box;
    std::shared_ptr<TemperatureController> m_controller;  

    Real m_t_air_interior = 0;
    Real m_t_sky4 = 0;
    Real m_t_air_exterior  = 0;
    Real m_t_air_exterior4 = 0;
    Real m_sigma = 5.6697E-8;  // Stefan-Boltzmann constant

};

class SolarThermalPostProcessor : public physics::PostProcessorBase
{
  public:
    explicit SolarThermalPostProcessor(std::shared_ptr<SourceTermSolarHeating> src_term) :
      m_src_term(src_term)
    {}

    // returns number of values this postprocessor returns
    int numValues() const override { return 2; }

    std::vector<std::string> getNames() const override { return {"solar_thermal_flux", "solar_thermal_controller"}; }

     std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override
     {
       return {m_src_term->getTotalFlux(), m_src_term->getControllerFactor()};
     }

  private:
    std::shared_ptr<SourceTermSolarHeating> m_src_term;
};

}

#endif