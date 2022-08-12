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

namespace Heat {

namespace impl {

struct RKStageState
{
  Real t         = 0;
  Real temp      = 0;
  Real dtemp_dt  = 0;
  Real flux      = 0;      // all fluxes except HVAC flux
  Real hvac_flux = 0; // HVAC flux
};

}

class InteriorAirTemperatureUpdator
{
  public:

    InteriorAirTemperatureUpdator(Real min_temp, Real max_temp, Real rho_cp, Real air_volume, 
                                  std::shared_ptr<AirLeakageModel> air_leakage, 
                                  std::shared_ptr<AirLeakageModel> ventilation,
                                  std::shared_ptr<InteriorLoads> interior_loads,
                                  std::shared_ptr<WindowConductionModel> window_model,
                                  Real initial_temp) :
      m_min_temp(min_temp),
      m_max_temp(max_temp),
      m_rho_cp(rho_cp),
      m_air_volume(air_volume),
      m_air_leakage(air_leakage),
      m_ventilation(ventilation),
      m_interior_loads(interior_loads),
      m_window_model(window_model),
      m_interior_temp(initial_temp),
      m_interior_temp_prev(initial_temp)
    {
      assertAlways(max_temp >= min_temp, "max temperature must be >= min temperature");
    }

    void initialize(HeatEquationSolar* heat_eqn_solar, const std::vector<NeumannBCPtr>& interior_bcs, DiscVectorPtr sol_vec, Real t_start)
    {
      m_t_prev = t_start;
      m_heat_eqn = heat_eqn_solar;
      m_bcs = interior_bcs;
      m_sol_vec_prev = makeDiscVector(heat_eqn_solar->getDiscretization());
      m_sol_vec_stage = makeDiscVector(heat_eqn_solar->getDiscretization());
      *m_sol_vec_prev = *sol_vec;
      computeInitialHVACFlux(sol_vec, t_start);
      m_net_flux_prev = computeNetFlux(sol_vec, t_start) + m_hvac_flux;
    }

    void updateTemperature(DiscVectorPtr sol_vec_np1, Real t)
    {
      assertAlways(t > m_t_prev, "new time must be > previous time");
      m_sol_vec_current = sol_vec_np1;
      m_t_current = t;
      computeRK();
      std::cout << "new temperature = " << m_interior_temp << std::endl;
      /*

      m_heat_eqn->setTimeParameters(t);
      Real flux_np1 = computeNetFlux(sol_vec_np1, t);

      Real delta_t = (t - m_t_prev) * 3600;
      Real fac = delta_t/(2 * m_rho_cp * m_air_volume);

      m_interior_temp = fac * (flux_np1 + m_net_flux_prev) + m_interior_temp_prev;
      std::cout << "new temperature = " << m_interior_temp << std::endl;

      if (m_interior_temp > m_max_temp)
        enforceTemperatureLimit(m_max_temp, flux_np1, delta_t);
      else if (m_interior_temp < m_min_temp)
        enforceTemperatureLimit(m_min_temp, flux_np1, delta_t);
      else
        m_hvac_flux = 0;

      m_net_flux_current = flux_np1 + m_hvac_flux;
      */
    }

    void startNewTimestep(DiscVectorPtr sol_vec_prev, Real t_prev)
    {
      std::cout << "starting new timestep" << std::endl;
      updateTemperature(sol_vec_prev, t_prev);

      m_net_flux_prev = m_net_flux_current;
      m_t_prev              = t_prev;
      m_interior_temp_prev  = m_interior_temp;
      *m_sol_vec_prev       = *sol_vec_prev;
    }

    Real getTemperature() const { return m_interior_temp; }

    Real getHVACFlux() const { return m_hvac_flux; }
    
  private:
    void computeRK()
    {
      Real delta_t = (m_t_current - m_t_prev) * 3600;

      auto state1 = computeRKStage(0, 0, m_interior_temp_prev);
      auto state2 = computeRKStage(0.5, 0.5, state1.dtemp_dt);
      auto state3 = computeRKStage(0.5, 0.5, state2.dtemp_dt);
      auto state4 = computeRKStage(1, 1, state3.dtemp_dt);

      m_interior_temp = m_interior_temp_prev + delta_t*(state1.dtemp_dt + 2 * state2.dtemp_dt + 2 * state3.dtemp_dt + state4.dtemp_dt)/6;
      m_net_flux_current = state4.flux + state4.hvac_flux;  //TODO: is this right, or should we recompute based in the new m_interior_temp;
    }

    impl::RKStageState computeRKStage(Real a_coeff, Real c_coeff, Real dtemp_dt_previous_stage)
    {
      impl::RKStageState state;
      Real delta_t_hours = m_t_current - m_t_prev;
      Real delta_t_seconds = delta_t_hours * 3600;

      state.t    = m_t_prev + c_coeff * delta_t_hours;
      state.temp = m_interior_temp_prev + a_coeff * delta_t_seconds * dtemp_dt_previous_stage;

      auto& u_prev_vec = m_sol_vec_prev->getVector();
      auto& u_curr_vec = m_sol_vec_current->getVector();
      auto& u_vec_stage = m_sol_vec_stage->getVector();
      for (int i=0; i < m_sol_vec_stage->getNumDofs(); ++i)
        u_vec_stage[i] = u_prev_vec[i] + c_coeff * (u_curr_vec[i] - u_prev_vec[i]);
      m_sol_vec_stage->markVectorModified();

      state.flux = computeNetFlux(m_sol_vec_stage, state.t);
      state.hvac_flux = 0; //TODO: need to enforce temperature limits

      state.dtemp_dt = state.flux/(m_rho_cp * m_air_volume);

      return state;
    }

    // computes HVAC flux such that the air temperature is the prescribed value
    void enforceTemperatureLimit(Real temp_limit, Real flux_np1, Real delta_t)
    {
      Real fac = delta_t/(2 * m_rho_cp * m_air_volume);
      m_hvac_flux = (temp_limit - m_interior_temp_prev) / fac - m_net_flux_prev - flux_np1;
      m_interior_temp = temp_limit;
    }

    // computes the initial HVAC flux such that the temperature is constant (ie the net flux is zero)
    void computeInitialHVACFlux(DiscVectorPtr sol_vec, Real t)
    {
      m_heat_eqn->setTimeParameters(t);
      m_hvac_flux = -computeNetFlux(sol_vec, t);
    }

    // flux from all sources except HVAC
    Real computeNetFlux(DiscVectorPtr sol_vec, Real t)
    {
      m_heat_eqn->setTimeParameters(t);

      Real flux = 0;
      for (auto& bc : m_bcs )
      {
        // BCs are defined such that flux into the wall is positive, so
        // the flux into the air is the negative
        flux -= computeFlux(bc, sol_vec, t);
      }

      flux -= m_air_leakage->computeAirLeakagePower(m_interior_temp, m_heat_eqn->getEnvironmentData().air_temp);
      flux -= m_ventilation->computeAirLeakagePower(m_interior_temp, m_heat_eqn->getEnvironmentData().air_temp);
      flux += m_interior_loads->computeLoadPower();
      flux += m_window_model->computeConductionPower(m_interior_temp, m_heat_eqn->getEnvironmentData().air_temp);

      return flux;
    }

    Real computeFlux(NeumannBCPtr bc, DiscVectorPtr sol_vec, Real t)
    {
      auto surf = bc->getSurfDisc();
      ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
      ArrayType<Real, 2> flux_vals(boost::extents[surf->getNumQuadPtsPerFace()][3]);
      std::vector<Real> flux_vals_vec(3*surf->getNumQuadPtsPerFace());

      if (!sol_vec->isArrayCurrent())
        sol_vec->syncVectorToArray();

      Real flux = 0;
      for (int face=0; face < surf->getNumFaces(); ++face)
      {
        auto& face_spec = surf->face_group.faces[face];
        auto& u_arr = sol_vec->getArray(face_spec.vol_group);
        auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];

        surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);

        bc->getValue(face, t, u_quad.data(), flux_vals_vec.data());
        for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
          for (int d=0; d < 3; ++d)
            flux_vals[i][d] = flux_vals_vec[d*surf->getNumQuadPtsPerFace() + i];

        flux += surf->getFaceWeight(face) * integrateFaceVector(surf, face, flux_vals);
      }

      Real global_flux = 0;
      MPI_Allreduce(&flux, &global_flux, 1, REAL_MPI_DATATYPE, MPI_SUM, MPI_COMM_WORLD);

      return global_flux;
    }

    // parameters
    Real m_min_temp;
    Real m_max_temp;
    Real m_rho_cp;
    Real m_air_volume;

    // sub-models
    HeatEquationSolar* m_heat_eqn;
    std::shared_ptr<AirLeakageModel> m_air_leakage;
    std::shared_ptr<AirLeakageModel> m_ventilation;
    std::shared_ptr<InteriorLoads> m_interior_loads;
    std::shared_ptr<WindowConductionModel> m_window_model;
    std::vector<NeumannBCPtr> m_bcs;

    // working state
    Real m_t_prev;
    Real m_t_current;
    DiscVectorPtr m_sol_vec_prev;     // solution at previous timestep
    DiscVectorPtr m_sol_vec_current;  // solution at current timestep
    DiscVectorPtr m_sol_vec_stage;    // vector used by RK method
    Real m_hvac_flux          = 0;
    Real m_interior_temp;
    Real m_interior_temp_prev;
    Real m_net_flux_prev      = 0;
    Real m_net_flux_current   = 0;
};

}

#endif