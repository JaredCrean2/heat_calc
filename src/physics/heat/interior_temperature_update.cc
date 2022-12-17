#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/HeatEquationSolar.h"

namespace Heat {

Real InteriorAirTemperatureUpdator::computeNetFlux(DiscVectorPtr sol_vec, Real interior_temp, Real t)
{
  Real load_flux = computeLoadFlux(sol_vec, interior_temp, t);

  if (interior_temp > m_max_temp)
    enforceTemperatureLimit(m_max_temp, interior_temp, load_flux);
  else if (interior_temp < m_min_temp)
    enforceTemperatureLimit(m_min_temp, interior_temp, load_flux);
  else
    m_hvac_flux = 0;

  std::cout << "load_flux = " << load_flux << ", hvac flux = " << m_hvac_flux << std::endl;

  return load_flux + m_hvac_flux;
}

Real InteriorAirTemperatureUpdator::computeNetFluxJacobian(DiscVectorPtr sol_vec, Real interior_temp, Real t)
{
  Real load_flux_dot;
  Real load_flux = computeLoadFluxDotTair(sol_vec, interior_temp, t, load_flux_dot);
  Real hvac_flux_dot;
  if (interior_temp > m_max_temp)
    hvac_flux_dot = enforceTemperatureLimitDotTair(m_max_temp, interior_temp, load_flux, load_flux_dot);
  else if (interior_temp < m_min_temp)
    hvac_flux_dot = enforceTemperatureLimitDotTair(m_min_temp, interior_temp, load_flux, load_flux_dot);
  else
    hvac_flux_dot = 0;  

  return load_flux_dot + hvac_flux_dot;
}

// compute d (computeNetFlux)/d (finite element solution)
void InteriorAirTemperatureUpdator::computeNetFlux_rev(DiscVectorPtr sol_vec, Real interior_temp, Real t, DiscVectorPtr sol_vec_bar, Real net_flux_bar)
{
  Real load_flux = computeLoadFlux(sol_vec, interior_temp, t);

  if (interior_temp > m_max_temp)
    enforceTemperatureLimit(m_max_temp, interior_temp, load_flux);
  else if (interior_temp < m_min_temp)
    enforceTemperatureLimit(m_min_temp, interior_temp, load_flux);
  else
    m_hvac_flux = 0;


  //-----------------------------------------------------------------------
  // Reverse mode
  //Real net_flux =  load_flux + m_hvac_flux;
  Real load_flux_bar = net_flux_bar;
  Real hvac_flux_bar = net_flux_bar;

  Real interior_temp_bar = 0;
  if (interior_temp > m_max_temp)
    enforceTemperatureLimit_rev(m_max_temp, interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
  else if (interior_temp < m_min_temp)
    enforceTemperatureLimit_rev(m_min_temp, interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
  
  computeLoadFlux_rev(sol_vec, sol_vec_bar, interior_temp, interior_temp_bar, t, load_flux_bar);
}


void InteriorAirTemperatureUpdator::enforceTemperatureLimit(Real temp_limit, Real interior_temp, Real load_flux)
{
  m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
}

Real InteriorAirTemperatureUpdator::enforceTemperatureLimitDotTair(Real temp_limit, Real interior_temp, Real load_flux, Real load_flux_dot)
{
  //m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  Real hvac_flux_dot = m_rho_cp * m_air_volume / m_hvac_restore_time - load_flux_dot;

  return hvac_flux_dot;
}

void InteriorAirTemperatureUpdator::enforceTemperatureLimit_rev(Real temp_limit, Real interior_temp, Real& interior_temp_bar, Real load_flux, 
                                  Real& load_flux_bar, Real hvac_flux_bar)
{
  //Real hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  interior_temp_bar += -m_rho_cp * m_air_volume * hvac_flux_bar / m_hvac_restore_time;
  load_flux_bar     -= hvac_flux_bar;
}


Real InteriorAirTemperatureUpdator::computeLoadFlux(DiscVectorPtr sol_vec, Real interior_temp, Real t)
{
  m_heat_eqn->setTimeParameters(t);

  Real flux = 0;
  for (auto& bc : m_bcs )
  {
    // BCs are defined such that flux into the wall is positive, so
    // the flux into the air is the negative
    flux -= computeBCFlux(bc, sol_vec, t);
  }

  flux -= m_air_leakage->computeAirLeakagePower(interior_temp, m_heat_eqn->getEnvironmentData().air_temp);
  flux -= m_ventilation->computeAirLeakagePower(interior_temp, m_heat_eqn->getEnvironmentData().air_temp);
  flux += m_interior_loads->computeLoadPower();
  flux += m_window_model->computeConductionPower(interior_temp, m_heat_eqn->getEnvironmentData().air_temp);

  return flux;
}

Real InteriorAirTemperatureUpdator::computeLoadFluxDotTair(DiscVectorPtr sol_vec, Real interior_temp, Real t, Real& flux_dot)
{
  flux_dot = 0;
  m_heat_eqn->setTimeParameters(t);

  Real flux = 0, flux_dot_tmp = 0;
  for (auto& bc : m_bcs )
  {
    // BCs are defined such that flux into the wall is positive, so
    // the flux into the air is the negative
    flux     -= computeBCFluxDotTair(bc, sol_vec, t, flux_dot_tmp);
    flux_dot -= flux_dot_tmp;
  }

  flux     -= m_air_leakage->computeAirLeakagePowerDot(interior_temp, m_heat_eqn->getEnvironmentData().air_temp, flux_dot_tmp);
  flux_dot -= flux_dot_tmp;

  flux     -= m_ventilation->computeAirLeakagePowerDot(interior_temp, m_heat_eqn->getEnvironmentData().air_temp, flux_dot_tmp);
  flux_dot -= flux_dot_tmp;

  flux += m_interior_loads->computeLoadPower();

  flux     += m_window_model->computeConductionPowerDot(interior_temp, m_heat_eqn->getEnvironmentData().air_temp, flux_dot_tmp);
  flux_dot += flux_dot_tmp;

  return flux;
}

void InteriorAirTemperatureUpdator::computeLoadFlux_rev(DiscVectorPtr sol_vec, DiscVectorPtr sol_vec_bar, Real interior_temp, Real& interior_temp_bar, Real t, Real flux_bar)
{
  m_heat_eqn->setTimeParameters(t);
/*
  Real flux = 0;
  for (auto& bc : m_bcs )
  {
    // BCs are defined such that flux into the wall is positive, so
    // the flux into the air is the negative
    flux -= computeBCFlux(bc, sol_vec, t);
  }

  flux -= m_air_leakage->computeAirLeakagePower(interior_temp, m_heat_eqn->getEnvironmentData().air_temp);
  flux -= m_ventilation->computeAirLeakagePower(interior_temp, m_heat_eqn->getEnvironmentData().air_temp);
  flux += m_interior_loads->computeLoadPower();
  flux += m_window_model->computeConductionPower(interior_temp, m_heat_eqn->getEnvironmentData().air_temp);
*/
  //-----------------------------------------------------------------------
  
  //TODO: interior_temp_bar is not used, don't compute it
  for (auto& bc : m_bcs )
  {
    // BCs are defined such that flux into the wall is positive, so
    // the flux into the air is the negative
    computeBCFlux_rev(bc, sol_vec, sol_vec_bar, t, -flux_bar);
  }

}


Real InteriorAirTemperatureUpdator::computeBCFlux(NeumannBCPtr bc, DiscVectorPtr sol_vec, Real t)
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
    {
      for (int d=0; d < 3; ++d)
        flux_vals[i][d] = flux_vals_vec[d*surf->getNumQuadPtsPerFace() + i];
    }

    flux += surf->getFaceWeight(face) * integrateFaceVector(surf, face, flux_vals);
  }

  Real global_flux = 0;
  MPI_Allreduce(&flux, &global_flux, 1, REAL_MPI_DATATYPE, MPI_SUM, MPI_COMM_WORLD);

  return global_flux;
}

Real InteriorAirTemperatureUpdator::computeBCFluxDotTair(NeumannBCPtr bc, DiscVectorPtr sol_vec, Real t, Real& flux_dot)
{
  auto bc_air_wind_sky = std::dynamic_pointer_cast<AirWindSkyNeumannBC>(bc);
  assertAlways(!!bc_air_wind_sky, "could not cast to AirWindSkyNeumannBC");

  auto surf = bc->getSurfDisc();
  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  ArrayType<Real, 2> flux_vals(boost::extents[surf->getNumQuadPtsPerFace()][3]);
  ArrayType<Real, 2> flux_vals_dot(boost::extents[surf->getNumQuadPtsPerFace()][3]);

  std::vector<Real> flux_vals_vec(3*surf->getNumQuadPtsPerFace());
  std::vector<Real> flux_vals_vec_dot(3*surf->getNumQuadPtsPerFace());

  if (!sol_vec->isArrayCurrent())
    sol_vec->syncVectorToArray();

  Real flux = 0;
  flux_dot = 0;
  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    auto& face_spec = surf->face_group.faces[face];
    auto& u_arr = sol_vec->getArray(face_spec.vol_group);
    auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];

    surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);

    bc_air_wind_sky->getValuedTair(face, t, u_quad.data(), flux_vals_vec.data(), flux_vals_vec_dot.data());
    for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
      for (int d=0; d < 3; ++d)
      {
        flux_vals[i][d] = flux_vals_vec[d*surf->getNumQuadPtsPerFace() + i];
        flux_vals_dot[i][d] = flux_vals_vec_dot[d*surf->getNumQuadPtsPerFace() + i];
      }

    flux += surf->getFaceWeight(face) * integrateFaceVector(surf, face, flux_vals);
    flux_dot += surf->getFaceWeight(face) * integrateFaceVector(surf, face, flux_vals_dot);
  }

  Real local_flux_and_dot[2]  = {flux, flux_dot};
  Real global_flux_and_dot[2] = {0, 0};
  MPI_Allreduce(local_flux_and_dot, global_flux_and_dot, 2, REAL_MPI_DATATYPE, MPI_SUM, MPI_COMM_WORLD);

  flux_dot = global_flux_and_dot[1];
  return global_flux_and_dot[0];
}


void InteriorAirTemperatureUpdator::computeBCFlux_rev(NeumannBCPtr bc, DiscVectorPtr sol_vec, DiscVectorPtr sol_vec_bar, Real t, Real global_flux_bar)
{
  auto bc_air_wind_sky = std::dynamic_pointer_cast<AirWindSkyNeumannBC>(bc);

  auto surf = bc->getSurfDisc();

  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  ArrayType<Real, 1> u_quad_bar(boost::extents[surf->getNumQuadPtsPerFace()]);

  ArrayType<Real, 2> flux_vals(boost::extents[surf->getNumQuadPtsPerFace()][3]);
  ArrayType<Real, 2> flux_vals_bar(boost::extents[surf->getNumQuadPtsPerFace()][3]);

  std::vector<Real> flux_vals_vec(3*surf->getNumQuadPtsPerFace());
  std::vector<Real> flux_vals_vec_bar(3*surf->getNumQuadPtsPerFace());


  if (!sol_vec->isArrayCurrent())
    sol_vec->syncVectorToArray();

  if (!sol_vec_bar->isArrayCurrent())
    sol_vec->syncVectorToArray();

  //-----------------------------------------------------------------------
  //TODO: zero out arrays
  Real flux_bar = global_flux_bar;
  //Real flux = 0;
  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    auto& face_spec = surf->face_group.faces[face];
    auto& u_arr = sol_vec->getArray(face_spec.vol_group);
    auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];

    auto& u_arr_bar = sol_vec_bar->getArray(face_spec.vol_group);
    auto u_el_bar = u_arr_bar[boost::indices[face_spec.el_group][range()]];

    surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);

    //---------------------------------------------------------------------
    zeroMatrix(flux_vals_bar);

    integrateFaceVector_rev(surf, face, flux_vals_bar, surf->getFaceWeight(face)*flux_bar);

    for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
      for (int d=0; d < 3; ++d)
        flux_vals_vec_bar[d*surf->getNumQuadPtsPerFace() + i] = flux_vals_bar[i][d]; 

    bc_air_wind_sky->getValue_rev(face, t, u_quad.data(), u_quad_bar.data(), flux_vals_vec_bar.data());

    surf->interp_vsq_flat[face_spec.face].interpolateVals_rev(u_el_bar, u_quad_bar);
  }

  sol_vec_bar->markArrayModified();

}

}