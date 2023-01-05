#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/HeatEquationSolar.h"
#include "spline.h"

namespace Heat {

void InteriorAirTemperatureUpdator::setExteriorTemperature(Real t_exterior)
{
  m_window_model->setExteriorTemperature(t_exterior);
  m_air_leakage->setExteriorTemperature(t_exterior);
  m_ventilation->setExteriorTemperature(t_exterior);
}


Real InteriorAirTemperatureUpdator::computeNetFlux(DiscVectorPtr sol_vec, Real interior_temp, Real t)
{
  Real load_flux = computeLoadFlux(sol_vec, interior_temp, t);
  m_hvac_flux = m_hvac_model->enforceTemperatureLimit(interior_temp, load_flux);

  std::cout << "load_flux = " << load_flux << ", hvac flux = " << m_hvac_flux << std::endl;

  return load_flux + m_hvac_flux;
}

Real InteriorAirTemperatureUpdator::computeNetFluxJacobian(DiscVectorPtr sol_vec, Real interior_temp, Real t)
{
  Real load_flux_dot;
  Real load_flux = computeLoadFluxDotTair(sol_vec, interior_temp, t, load_flux_dot);
  Real hvac_flux_dot = m_hvac_model->enforceTemperatureLimit_dot(interior_temp, 1, load_flux, load_flux_dot); 

  std::cout << "load_flux = " << load_flux << std::endl;
  std::cout << "load_flux_dot = " << load_flux_dot << ", hvac_flux_dot = " << hvac_flux_dot << std::endl;
  std::cout << "net flux dot = " << load_flux_dot + hvac_flux_dot << std::endl;

  return load_flux_dot + hvac_flux_dot;
}

// compute d (computeNetFlux)/d (finite element solution)
void InteriorAirTemperatureUpdator::computeNetFlux_rev(DiscVectorPtr sol_vec, Real interior_temp, Real t, DiscVectorPtr sol_vec_bar, Real net_flux_bar)
{
  Real load_flux = computeLoadFlux(sol_vec, interior_temp, t);
  m_hvac_flux = m_hvac_model->enforceTemperatureLimit(interior_temp, load_flux);

  //-----------------------------------------------------------------------
  // Reverse mode
  //Real net_flux =  load_flux + m_hvac_flux;
  Real load_flux_bar = net_flux_bar;
  Real hvac_flux_bar = net_flux_bar;

  Real interior_temp_bar = 0;
  m_hvac_model->enforceTemperatureLimit_rev(interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);  
  computeLoadFlux_rev(sol_vec, sol_vec_bar, interior_temp, interior_temp_bar, t, load_flux_bar);
}

/*
Real InteriorAirTemperatureUpdator::enforceTemperatureLimitStraightLine(Real temp_limit, Real interior_temp, Real load_flux)
{
  return m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
}

Real InteriorAirTemperatureUpdator::enforceTemperatureLimitStraightLineDotTair(Real temp_limit, Real interior_temp,
                                                                               Real load_flux, Real load_flux_dot)
{
  //m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  Real hvac_flux_dot = m_rho_cp * m_air_volume / m_hvac_restore_time - load_flux_dot;

  return hvac_flux_dot;
}

void InteriorAirTemperatureUpdator::enforceTemperatureLimitStraightLineDotTair_rev(Real temp_limit, Real interior_temp,
                                                                                   Real load_flux,
                                                                                   Real load_flux_dot, Real& load_flux_dot_bar,
                                                                                   Real hvac_flux_dot_bar)
{
  //m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  //Real hvac_flux_dot = m_rho_cp * m_air_volume / m_hvac_restore_time - load_flux_dot;

  load_flux_dot_bar -= hvac_flux_dot_bar;
}

void InteriorAirTemperatureUpdator::enforceTemperatureLimitStraightLine_rev(Real temp_limit, Real interior_temp,
                                                                            Real& interior_temp_bar, Real load_flux, 
                                                                            Real& load_flux_bar, Real hvac_flux_bar)
{
  //Real hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  interior_temp_bar += -m_rho_cp * m_air_volume * hvac_flux_bar / m_hvac_restore_time;
  load_flux_bar     -= hvac_flux_bar;
}

// setup spline to smooth the activation of the HVAC system, which causes
// convergence problems for Newtons method
std::array<Real, 6> InteriorAirTemperatureUpdator::getSplineParams(Real interior_temp, Real load_flux)
{
  Real expansion_factor = 0.1;
  Real delta_t = (m_max_temp - m_min_temp);
  Real t_upper = m_max_temp + expansion_factor * delta_t;
  Real t_lower = m_min_temp - expansion_factor * delta_t;
  Real load_flux_dot = 0; // We are computing the spline so the HVAC flux is a continuous function of the air temperature,
                          // so don't consider the derivative wrt load flux

  Real flux_upper     = enforceTemperatureLimitStraightLine(m_max_temp, t_upper, load_flux);
  Real flux_upper_dot = enforceTemperatureLimitStraightLineDotTair(m_max_temp, t_upper, load_flux, load_flux_dot);

  Real flux_lower     = enforceTemperatureLimitStraightLine(m_min_temp, t_lower, load_flux);
  Real flux_lower_dot = enforceTemperatureLimitStraightLineDotTair(m_min_temp, t_lower, load_flux, load_flux_dot);

  return {t_lower, t_upper, flux_lower, flux_lower_dot, flux_upper, flux_upper_dot};
}

void InteriorAirTemperatureUpdator::getSplineParams_rev(Real interior_temp, Real& interior_temp_bar,
                                                        Real load_flux, Real& load_flux_bar,
                                                        const std::array<Real, 6> params_bar)
{
  Real expansion_factor = 0.1;
  Real delta_t = (m_max_temp - m_min_temp);
  Real t_upper = m_max_temp + expansion_factor * delta_t;
  Real t_lower = m_min_temp - expansion_factor * delta_t;
  Real load_flux_dot = 0; // We are computing the spline so the HVAC flux is a continuous function of the air temperature,
                          // so don't consider the derivative wrt load flux

  //Real flux_upper     = enforceTemperatureLimitStraightLine(m_max_temp, t_upper, load_flux);
  //Real flux_upper_dot = enforceTemperatureLimitStraightLineDotTair(m_max_temp, t_upper, load_flux, load_flux_dot);

  //Real flux_lower     = enforceTemperatureLimitStraightLine(m_min_temp, t_lower, load_flux);
  //Real flux_lower_dot = enforceTemperatureLimitStraightLineDotTair(m_min_temp, t_lower, load_flux, load_flux_dot);

  Real load_flux_dot_bar = 0;
  enforceTemperatureLimitStraightLineDotTair_rev(m_min_temp, t_lower, load_flux, load_flux_dot, load_flux_dot_bar, params_bar[3]);
  enforceTemperatureLimitStraightLine_rev(m_min_temp, t_lower, interior_temp_bar, load_flux, load_flux_bar, params_bar[2]);

  enforceTemperatureLimitStraightLineDotTair_rev(m_max_temp, t_upper, load_flux, load_flux_dot, load_flux_dot_bar, params_bar[5]);
  enforceTemperatureLimitStraightLine_rev(m_max_temp, t_upper, interior_temp_bar, load_flux, load_flux_bar, params_bar[4]);
}

Real InteriorAirTemperatureUpdator::enforceTemperatureLimit(Real interior_temp, Real load_flux)
{
  std::array<Real, 6> params = getSplineParams(interior_temp, load_flux);

  if (interior_temp > params[1])
    return enforceTemperatureLimitStraightLine(m_max_temp, interior_temp, load_flux);
  else if (interior_temp < params[0])
    return enforceTemperatureLimitStraightLine(m_min_temp, interior_temp, load_flux);
  else
  {
    SingleCubicSpline spline;
    spline.setupSpline(params[0], params[1], params[2], params[3], params[4], params[5]);
    return spline.eval(interior_temp);
  }
}

Real InteriorAirTemperatureUpdator::enforceTemperatureLimitDotTair(Real interior_temp, Real load_flux, Real load_flux_dot)
{
  std::array<Real, 6> params = getSplineParams(interior_temp, load_flux);

  if (interior_temp > params[1])
    return enforceTemperatureLimitStraightLineDotTair(m_max_temp, interior_temp, load_flux, load_flux_dot);
  else if (interior_temp < params[0])
    return enforceTemperatureLimitStraightLineDotTair(m_min_temp, interior_temp, load_flux, load_flux_dot);
  else
  {
    SingleCubicSpline spline;
    spline.setupSpline(params[0], params[1], params[2], params[3], params[4], params[5]);
    Real hvac_flux_dot;
    spline.evalDot(interior_temp, 1, hvac_flux_dot);
    return hvac_flux_dot;
  }
}

void InteriorAirTemperatureUpdator::enforceTemperatureLimit_rev(Real interior_temp, Real& interior_temp_bar,
                                                                Real load_flux, Real& load_flux_bar,
                                                                Real hvac_flux_bar)
{
  std::array<Real, 6> params = getSplineParams(interior_temp, load_flux);

  if (interior_temp > params[1])
    enforceTemperatureLimitStraightLine_rev(m_max_temp, interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
  else if (interior_temp < params[0])
    enforceTemperatureLimitStraightLine_rev(m_min_temp, interior_temp, interior_temp_bar, load_flux,load_flux_bar, hvac_flux_bar);
  else
  {
    SingleCubicSpline spline;
    spline.setupSpline(params[0], params[1], params[2], params[3], params[4], params[5]);

    std::array<Real, 4> coeffs_bar;
    spline.evalRev(interior_temp, hvac_flux_bar, interior_temp_bar, coeffs_bar);

    std::array<Real, 6> params_bar = {0};
    spline.setupSplineRev(params[0], params[1], coeffs_bar, params_bar[2], params_bar[3], params_bar[4], params_bar[5]);

    getSplineParams_rev(interior_temp, interior_temp_bar, load_flux, load_flux_bar, params_bar);
  }
}

*/



Real InteriorAirTemperatureUpdator::computeLoadFlux(DiscVectorPtr sol_vec, Real interior_temp, Real t)
{
  m_heat_eqn->setTimeParameters(t, interior_temp);

  Real bc_flux = 0;
  /*
  std::cout << "number of bcs = " << m_bcs.size() << std::endl;
  for (auto& bc : m_bcs )
  {
    // BCs are defined such that flux into the wall is positive, so
    // the flux into the air is the negative
    bc_flux -= computeBCFlux(bc, sol_vec, t);  
  }
  */
/*
  Real air_leakage_flux   = -m_air_leakage->computeAirLeakagePower(interior_temp);
  Real ventilation_flux   = -m_ventilation->computeAirLeakagePower(interior_temp);
  Real interior_load_flux = m_interior_loads->computeLoadPower();
  Real window_flux        = m_window_model->computeConductionPower(interior_temp);


  std::cout << "bc_flux            = " << bc_flux << std::endl;
  std::cout << "air_leakage_flux   = " << air_leakage_flux << std::endl;
  std::cout << "ventilation_flux   = " << ventilation_flux << std::endl;
  std::cout << "interior_load_flux = " << interior_load_flux << std::endl;
  std::cout << "window_flux        = " << window_flux << std::endl;
*/
  return bc_flux; // + air_leakage_flux + ventilation_flux + interior_load_flux + window_flux;
}

Real InteriorAirTemperatureUpdator::computeLoadFluxDotTair(DiscVectorPtr sol_vec, Real interior_temp, Real t, Real& flux_dot)
{
  flux_dot = 0;
  m_heat_eqn->setTimeParameters(t, interior_temp);

  Real flux = 0;
  /*
  for (auto& bc : m_bcs )
  {

    // BCs are defined such that flux into the wall is positive, so
    // the flux into the air is the negative
    Real flux_dot_tmp = 0;
    flux     -= computeBCFluxDotTair(bc, sol_vec, t, flux_dot_tmp);
    flux_dot -= flux_dot_tmp;  
  }
  */
/*
  Real flux_dot_tmp = 0;
  flux     -= m_air_leakage->computeAirLeakagePowerDot(interior_temp, flux_dot_tmp);
  flux_dot -= flux_dot_tmp;

  flux_dot_tmp = 0;
  flux     -= m_ventilation->computeAirLeakagePowerDot(interior_temp, flux_dot_tmp);
  flux_dot -= flux_dot_tmp;

  flux += m_interior_loads->computeLoadPower();

  flux_dot_tmp = 0;
  flux     += m_window_model->computeConductionPowerDot(interior_temp, flux_dot_tmp);
  flux_dot += flux_dot_tmp;
*/
  return flux;
}

void InteriorAirTemperatureUpdator::computeLoadFlux_rev(DiscVectorPtr sol_vec, DiscVectorPtr sol_vec_bar, Real interior_temp,
                                                        Real& interior_temp_bar, Real t, Real flux_bar)
{
  m_heat_eqn->setTimeParameters(t, interior_temp);
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
  /*
  for (auto& bc : m_bcs )
  {
    // BCs are defined such that flux into the wall is positive, so
    // the flux into the air is the negative
    computeBCFlux_rev(bc, sol_vec, sol_vec_bar, t, -flux_bar);
  }
  */

  //TODO: what about other models? interior_temp_bar?
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
    for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
      u_quad_bar[i] = 0;

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