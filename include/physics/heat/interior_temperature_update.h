#ifndef PHYSICS_HEAT_INTERIOR_TEMPERATURE_UPDATE_H
#define PHYSICS_HEAT_INTERIOR_TEMPERATURE_UPDATE_H

#include "ProjectDefs.h"
#include "discretization/NeumannBC.h"
#include "discretization/disc_vector.h"
#include "HeatEquation.h"
#include "discretization/surface_discretization.h"
#include "air_leakage.h"
#include "interior_loads.h"
#include "window_conduction_model.h"
#include "physics/heat/helper_funcs.h"

namespace Heat {

class HeatEquationSolar;

class InteriorAirTemperatureUpdator
{
  public:

    InteriorAirTemperatureUpdator(Real min_temp, Real max_temp, Real rho_cp, Real air_volume, 
                                  std::shared_ptr<AirLeakageModel> air_leakage, 
                                  std::shared_ptr<AirLeakageModel> ventilation,
                                  std::shared_ptr<InteriorLoads> interior_loads,
                                  std::shared_ptr<WindowConductionModel> window_model,
                                  Real hvac_restore_time) :
      m_min_temp(min_temp),
      m_max_temp(max_temp),
      m_rho_cp(rho_cp),
      m_air_volume(air_volume),
      m_hvac_restore_time(hvac_restore_time),
      m_air_leakage(air_leakage),
      m_ventilation(ventilation),
      m_interior_loads(interior_loads),
      m_window_model(window_model)
    {
      assertAlways(max_temp >= min_temp, "max temperature must be >= min temperature");
    }

    void initialize(HeatEquationSolar* heat_eqn_solar, const std::vector<NeumannBCPtr>& interior_bcs)
    {
      m_heat_eqn = heat_eqn_solar;
      m_bcs = interior_bcs;
    }

    Real computeNetFlux(DiscVectorPtr sol_vec, Real interior_temp, Real t);

    // compute d (computeNetFlux)/d (interior_temp)
    Real computeNetFluxJacobian(DiscVectorPtr sol_vec, Real interior_temp, Real t);

    // compute d (computeNetFlux)/d (finite element solution)
    void computeNetFlux_rev(DiscVectorPtr sol_vec, Real interior_temp, Real t, DiscVectorPtr sol_vec_bar, Real net_flux_bar);


    Real getThermalMass() const { return m_rho_cp * m_air_volume; }

    Real getHVACFlux() const { return m_hvac_flux; }
    
  private:

    // computes HVAC flux such that the air temperature returns to the temp_limit within
    // m_hvac_restore_time (approximately)
    void enforceTemperatureLimit(Real temp_limit, Real interior_temp, Real load_flux);

    Real enforceTemperatureLimitDotTair(Real temp_limit, Real interior_temp, Real load_flux, Real load_flux_dot);

    void enforceTemperatureLimit_rev(Real temp_limit, Real interior_temp, Real& interior_temp_bar, Real load_flux, 
                                     Real& load_flux_bar, Real hvac_flux_bar);


    // flux from all sources except HVAC
    Real computeLoadFlux(DiscVectorPtr sol_vec, Real interior_temp, Real t);

    Real computeLoadFluxDotTair(DiscVectorPtr sol_vec, Real interior_temp, Real t, Real& flux_dot);

    void computeLoadFlux_rev(DiscVectorPtr sol_vec, DiscVectorPtr sol_vec_bar, Real interior_temp, Real& interior_temp_bar, Real t, Real flux_bar);


    Real computeBCFlux(NeumannBCPtr bc, DiscVectorPtr sol_vec, Real t);

    Real computeBCFluxDotTair(NeumannBCPtr bc, DiscVectorPtr sol_vec, Real t, Real& flux_dot);

    void computeBCFlux_rev(NeumannBCPtr bc, DiscVectorPtr sol_vec, DiscVectorPtr sol_vec_bar, Real t, Real global_flux_bar);

    // parameters
    Real m_min_temp;
    Real m_max_temp;
    Real m_rho_cp;
    Real m_air_volume;
    Real m_hvac_restore_time;  // after the temperature exceeds the bounds, the HVAC
                               // system will restore it within (approximately) this
                               // much time, in seconds

    // sub-models
    HeatEquationSolar* m_heat_eqn;
    std::shared_ptr<AirLeakageModel> m_air_leakage;
    std::shared_ptr<AirLeakageModel> m_ventilation;
    std::shared_ptr<InteriorLoads> m_interior_loads;
    std::shared_ptr<WindowConductionModel> m_window_model;
    std::vector<NeumannBCPtr> m_bcs;

    // working state
    DiscVectorPtr m_sol_vec_prev;     // solution at previous timestep
    DiscVectorPtr m_sol_vec_current;  // solution at current timestep
    DiscVectorPtr m_sol_vec_stage;    // vector used by RK method
    Real m_hvac_flux          = 0;
    //Real m_interior_temp;
    //Real m_interior_temp_prev;
    //Real m_net_flux_prev      = 0;
};

}

#endif