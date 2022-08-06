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
      m_heat_eqn = heat_eqn_solar;
      m_bcs = interior_bcs;
      m_sol_vec_prev = makeDiscVector(heat_eqn_solar->getDiscretization());
      computeInitialHVACFlux(sol_vec, t_start);
      startNewTimestep(sol_vec, t_start);
    }

    void updateTemperature(DiscVectorPtr sol_vec_np1, Real t)
    {
      assertAlways(t > m_t_prev, "new time must be > previous time");

      m_heat_eqn->setTimeParameters(t);
      Real flux_np1 = computeNetFlux(sol_vec_np1, t);

      Real delta_t = (t - m_t_prev) * 3600;
      std::cout << "flux_np1 total = " << flux_np1 * delta_t << std::endl;
      Real fac = delta_t/(2 * m_rho_cp * m_air_volume);

      m_interior_temp = fac * (flux_np1 + m_net_flux_prev) + m_interior_temp_prev;
      //std::cout << "m_interior_temp = " << m_interior_temp << std::endl;

      if (m_interior_temp > m_max_temp)
        enforceTemperatureLimit(m_max_temp, flux_np1, delta_t);
      else if (m_interior_temp < m_min_temp)
        enforceTemperatureLimit(m_min_temp, flux_np1, delta_t);
      else
        m_hvac_flux = 0;
    }

    void startNewTimestep(DiscVectorPtr sol_vec_prev, Real t_prev)
    {
      m_t_prev              = t_prev;
      m_interior_temp_prev  = m_interior_temp;
      *m_sol_vec_prev       = *sol_vec_prev;

      m_heat_eqn->setTimeParameters(t_prev);
      //TODO: is this the best way of doing things.  Maybe call updateTemperature instead?
      m_net_flux_prev = computeNetFlux(m_sol_vec_prev, t_prev) + m_hvac_flux;
    }

    Real getTemperature() const { return m_interior_temp; }

    Real getHVACFlux() const { return m_hvac_flux; }
    
  private:

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
    DiscVectorPtr m_sol_vec_prev;
    Real m_hvac_flux          = 0;
    Real m_interior_temp;
    Real m_interior_temp_prev;
    Real m_net_flux_prev      = 0;
};

}

#endif