#include "physics/heat/bc_defs.h"

namespace Heat {

//-----------------------------------------------------------------------------
// NewtonCooling

void NewtonCooling::getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
{
  auto& normals = m_surf->normals;
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    Real nx = normals[face][i][0];
    Real ny = normals[face][i][1];
    Real nz = normals[face][i][2];
    Real mag     = std::sqrt(nx*nx + ny*ny + nz*nz);
    Real delta_t = m_temp - sol_vals[i];
    Real term    = m_heat_transfer_coeff * delta_t / mag;

    flux_vals[i]                                      = nx * term;
    flux_vals[i + m_surf->getNumQuadPtsPerFace()]     = ny * term;
    flux_vals[i + 2 * m_surf->getNumQuadPtsPerFace()] = nz * term;
  }
}


void NewtonCooling::getValueDeriv(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals_deriv)
{
  auto& normals = m_surf->normals;
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    Real nx = normals[face][i][0];
    Real ny = normals[face][i][1];
    Real nz = normals[face][i][2];
    Real mag         = std::sqrt(nx*nx + ny*ny + nz*nz);
    Real delta_t_dot = -1;
    Real term        = m_heat_transfer_coeff * delta_t_dot / mag;

    flux_vals_deriv[i]                                      = nx * term;
    flux_vals_deriv[i + m_surf->getNumQuadPtsPerFace()]     = ny * term;
    flux_vals_deriv[i + 2 * m_surf->getNumQuadPtsPerFace()] = nz * term;
  }
}

//-----------------------------------------------------------------------------
void NewtonCoolingFromAir::getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
{
  auto& normals = m_surf->normals;
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    Real nx = normals[face][i][0];
    Real ny = normals[face][i][1];
    Real nz = normals[face][i][2];
    Real mag     = std::sqrt(nx*nx + ny*ny + nz*nz);
    Real delta_t = m_temp - sol_vals[i];
    Real term    = m_heat_transfer_coeff * delta_t / mag;

    flux_vals[i]                                      = nx * term;
    flux_vals[i + m_surf->getNumQuadPtsPerFace()]     = ny * term;
    flux_vals[i + 2 * m_surf->getNumQuadPtsPerFace()] = nz * term;
  }
}


void NewtonCoolingFromAir::getValueDeriv(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals_deriv)
{
  auto& normals = m_surf->normals;
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    Real nx = normals[face][i][0];
    Real ny = normals[face][i][1];
    Real nz = normals[face][i][2];
    Real mag         = std::sqrt(nx*nx + ny*ny + nz*nz);
    Real delta_t_dot = -1;
    Real term        = m_heat_transfer_coeff * delta_t_dot / mag;

    flux_vals_deriv[i]                                      = nx * term;
    flux_vals_deriv[i + m_surf->getNumQuadPtsPerFace()]     = ny * term;
    flux_vals_deriv[i + 2 * m_surf->getNumQuadPtsPerFace()] = nz * term;
  }
}

void NewtonCoolingFromAir::getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv)
{
  getValue(face, t, sol_vals, flux_vals);
  getValueDeriv(face, t, sol_vals, flux_vals_deriv);
  for (int i=0; i < 3*m_surf->getNumQuadPtsPerFace(); ++i)
    flux_vals_deriv[i] *= -1;
}

void NewtonCoolingFromAir::getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar)
{
  auto& normals = m_surf->normals;
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    Real nx = normals[face][i][0];
    Real ny = normals[face][i][1];
    Real nz = normals[face][i][2];
    Real mag     = std::sqrt(nx*nx + ny*ny + nz*nz);
    //Real delta_t = m_temp - sol_vals[i];
    //Real term    = m_heat_transfer_coeff * delta_t / mag;

    //flux_vals[i]                                      = nx * term;
    //flux_vals[i + m_surf->getNumQuadPtsPerFace()]     = ny * term;
    //flux_vals[i + 2 * m_surf->getNumQuadPtsPerFace()] = nz * term;

    Real term_bar = 0;
    term_bar += nx * flux_vals_bar[i];
    term_bar += ny * flux_vals_bar[i + m_surf->getNumQuadPtsPerFace()];
    term_bar += nz * flux_vals_bar[i + 2 * m_surf->getNumQuadPtsPerFace()];

    Real delta_t_bar = m_heat_transfer_coeff * term_bar / mag;
    sol_vals_bar[i] = - delta_t_bar;
  }
}

//-----------------------------------------------------------------------------
// TarpBC

void TarpBC::getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
{
  getQuadNodeCoords(face, m_quad_coords);
  Real air_temp = m_tarp.getAirTemperature();
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> pt{m_quad_coords[i][0], m_quad_coords[i][1], m_quad_coords[i][2]};
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));
    Real wall_temp = sol_vals[i];

    Real heat_transfer_coeff = m_tarp.computeHeatTransferCoeff(wall_temp, pt, unit_normal);

    for (int d=0; d < 3; ++d)
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * heat_transfer_coeff * (air_temp - wall_temp);
  }
}

void TarpBC::getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) 
{
  getQuadNodeCoords(face, m_quad_coords);
  Real air_temp = m_tarp.getAirTemperature();
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> pt{m_quad_coords[i][0], m_quad_coords[i][1], m_quad_coords[i][2]};
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));
    Real wall_temp = sol_vals[i];

    Real heat_transfer_coeff_dot;
    Real heat_transfer_coeff = m_tarp.computeHeatTransferCoeffdTwall(wall_temp, pt, unit_normal, heat_transfer_coeff_dot);

    for (int d=0; d < 3; ++d)
    {
      //flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * heat_transfer_coeff * (wall_temp - air_temp);
      flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * heat_transfer_coeff_dot * (air_temp - wall_temp) -
                                                                unit_normal[d] * heat_transfer_coeff;
    }
  }
}

void TarpBC::getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv)
{
  getQuadNodeCoords(face, m_quad_coords);
  Real air_temp = m_tarp.getAirTemperature();
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> pt{m_quad_coords[i][0], m_quad_coords[i][1], m_quad_coords[i][2]};
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));
    Real wall_temp = sol_vals[i];

    Real heat_transfer_coeff_dot;
    Real heat_transfer_coeff = m_tarp.computeHeatTransferCoeffdTair(wall_temp, pt, unit_normal, heat_transfer_coeff_dot);

    for (int d=0; d < 3; ++d)
    {
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i]       = unit_normal[d] * heat_transfer_coeff * (air_temp - wall_temp);
      flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * heat_transfer_coeff_dot * (air_temp - wall_temp) +
                                                                unit_normal[d] * heat_transfer_coeff;
    }
  }
}

void TarpBC::getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar)
{
  getQuadNodeCoords(face, m_quad_coords);
  Real air_temp = m_tarp.getAirTemperature();
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> pt{m_quad_coords[i][0], m_quad_coords[i][1], m_quad_coords[i][2]};
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));
    Real wall_temp = sol_vals[i];
    Real wall_temp_bar = 0;

    Real heat_transfer_coeff_bar = 0, heat_transfer_coeff_dot;
    Real heat_transfer_coeff = m_tarp.computeHeatTransferCoeffdTwall(wall_temp, pt, unit_normal, heat_transfer_coeff_dot);

    for (int d=0; d < 3; ++d)
    {
      //flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * heat_transfer_coeff * (air_temp - wall_temp);
      Real flux_val_bar = flux_vals_bar[d * m_surf->getNumQuadPtsPerFace() + i];
      wall_temp_bar += -unit_normal[d] * heat_transfer_coeff * flux_val_bar;
      heat_transfer_coeff_bar += unit_normal[d] * (air_temp - wall_temp) * flux_val_bar;
    }

    wall_temp_bar += heat_transfer_coeff_dot * heat_transfer_coeff_bar;
    sol_vals_bar[i] = wall_temp_bar;
  }
}

//-----------------------------------------------------------------------------
// SkyRadiationBC

void SkyRadiationBC::getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    auto unit_normal = getUnitNormal(face, i);
    Real flux = m_model.computeFlux(sol_vals[i], unit_normal);
    for (int d=0; d < 3; ++d)
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
  }
}

void SkyRadiationBC::getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    auto unit_normal = getUnitNormal(face, i);
    Real flux_dot = 0;
    m_model.computeFluxdTwall(sol_vals[i], unit_normal, flux_dot);
    for (int d=0; d < 3; ++d)
      flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux_dot;
  }
}

void SkyRadiationBC::getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    auto unit_normal = getUnitNormal(face, i);
    Real flux_dot = 0;
    Real flux = m_model.computeFluxdTair(sol_vals[i], unit_normal, flux_dot);
    for (int d=0; d < 3; ++d)
    {
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
      flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux_dot;
    }
  }
}

void SkyRadiationBC::getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    auto unit_normal = getUnitNormal(face, i);
    Real flux_dot;
    m_model.computeFluxdTwall(sol_vals[i], unit_normal, flux_dot);
    Real flux_bar = 0;
    for (int d=0; d < 3; ++d)
    {
      //flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
      flux_bar += unit_normal[d] * flux_vals_bar[d * m_surf->getNumQuadPtsPerFace() + i];
    }

    sol_vals_bar[i] = flux_dot * flux_bar;
  }
}

//-----------------------------------------------------------------------------
// SolarRadiationBC

void SolarRadiationBC::getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    auto unit_normal = getUnitNormal(face, i);
    Real flux = m_model.computeFlux(unit_normal);
    if (flux < 0)
      throw std::runtime_error("flux cannot be negative");
    for (int d=0; d < 3; ++d)
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
  }
}

void SolarRadiationBC::getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv)
{
  for (int i=0; i < 3*m_surf->getNumQuadPtsPerFace(); ++i)
    flux_vals_deriv[i] = 0;
}

void SolarRadiationBC::getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
    sol_vals_bar[i] = 0;
}

//-----------------------------------------------------------------------------
// SimpleConvectionBC
void SimpleConvectionBC::getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));

    for (int d=0; d < 3; ++d)
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * m_heat_transfer_coeff * (m_air_temp - sol_vals[i]);
  }
}

void SimpleConvectionBC::getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));

    for (int d=0; d < 3; ++d)
      flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = -unit_normal[d] * m_heat_transfer_coeff;
  }
}

void SimpleConvectionBC::getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));

    for (int d=0; d < 3; ++d)
    {
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * m_heat_transfer_coeff * (m_air_temp - sol_vals[i]);
      flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * m_heat_transfer_coeff;
    }
  }
}

void SimpleConvectionBC::getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar)
{
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));
    sol_vals_bar[i] = 0;

    for (int d=0; d < 3; ++d)
    {
      //flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * m_heat_transfer_coeff * (m_air_temp - sol_vals[i]);
      sol_vals_bar[i] += -unit_normal[d] * m_heat_transfer_coeff * flux_vals_bar[d * m_surf->getNumQuadPtsPerFace() + i];
    }
  }
}

//-----------------------------------------------------------------------------
// FloorRadiationBC

void FloorRadiationBC::getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
{
  Real flux = m_model.computeFlux();
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));

    for (int d=0; d < 3; ++d)
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
  }
}


void FloorRadiationBC::getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv)
{
  Real flux = m_model.computeFlux();
  for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
  {
    std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
    auto unit_normal = normal / std::sqrt(dot(normal, normal));

    for (int d=0; d < 3; ++d)
    {
      flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
      flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = 0;
    }
  }
}


void FloorRadiationBC::getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar)
{
  // nothing to do here
}


//-----------------------------------------------------------------------------
// CombinedAirWindSkyBCs

namespace
{

bool anyBcsNonlinear(const std::vector<std::shared_ptr<AirWindSkyNeumannBC>>& bcs)
{
  bool is_nonlinear = false;
  for (auto& bc : bcs)
    is_nonlinear = is_nonlinear || bc->isNonlinear();

  return is_nonlinear;
}

}

CombinedAirWindSkyNeumannBC::CombinedAirWindSkyNeumannBC(std::vector<std::shared_ptr<AirWindSkyNeumannBC>> bcs) :
  AirWindSkyNeumannBC(bcs[0]->getSurfDisc(), anyBcsNonlinear(bcs), "combined"),
  m_bcs(bcs)
{
  auto surf0 = bcs[0]->getSurfDisc();
  for (size_t i=1; i < bcs.size(); ++i)
    if (surf0 != bcs[i]->getSurfDisc())
      throw std::runtime_error("BCs have different surfaces");
}



void CombinedAirWindSkyNeumannBC::setAirTemperature(Real temp)
{
  for (auto& bc : m_bcs)
    bc->setAirTemperature(temp);
}

Real CombinedAirWindSkyNeumannBC::getAirTemperature() const
{ 
  return m_bcs[0]->getAirTemperature();
}

void CombinedAirWindSkyNeumannBC::setAirSpeed(Real velocity)
{
  for (auto& bc : m_bcs)
    bc->setAirSpeed(velocity);
}

void CombinedAirWindSkyNeumannBC::setAirDirection(std::array<Real, 3> direction)
{
  for (auto& bc : m_bcs)
    bc->setAirDirection(direction);
}

void CombinedAirWindSkyNeumannBC::setIRHorizontalRadiation(Real flux)
{
  for (auto& bc : m_bcs)
    bc->setIRHorizontalRadiation(flux);
}

void CombinedAirWindSkyNeumannBC::setDirectNormalRadiation(Real flux)
{
  for (auto& bc : m_bcs)
    bc->setDirectNormalRadiation(flux);
}

void CombinedAirWindSkyNeumannBC::setDiffuseRadiation(Real flux)
{
  for (auto& bc : m_bcs)
    bc->setDiffuseRadiation(flux);
}    

void CombinedAirWindSkyNeumannBC::setSolarDirection(const DirectionCosines& cosines)
{
  for (auto& bc : m_bcs)
    bc->setSolarDirection(cosines);
}

void CombinedAirWindSkyNeumannBC::getValue(const Index face, const Real t, const Real* sol_vals, Real* flux_vals)
{
  int npts = getSurfDisc()->getNumQuadPtsPerFace()*3;
  std::vector<Real> flux_vals_tmp(npts, 0);

  for (int i=0; i < npts; ++i)
    flux_vals[i] = 0;

  //int i=0;
  for (auto& bc : m_bcs)
  {
    bc->getValue(face, t, sol_vals, flux_vals_tmp.data());
    updateAndZero(flux_vals, flux_vals_tmp.data(), npts);
    //++i;
  }
}

void CombinedAirWindSkyNeumannBC::getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv)
{
  int npts = getSurfDisc()->getNumQuadPtsPerFace()*3;
  std::vector<Real> flux_vals_deriv_tmp(npts, 0);

  for (int i=0; i < npts; ++i)
    flux_vals_deriv[i] = 0;

  for (auto& bc : m_bcs)
  {
    bc->getValueDeriv(face, t, sol_vals, flux_vals_deriv_tmp.data());
    updateAndZero(flux_vals_deriv, flux_vals_deriv_tmp.data(), npts);
  }
}

// compute derivative of flux_vals wrt air temperature
void CombinedAirWindSkyNeumannBC::getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv)
{
  int npts = getSurfDisc()->getNumQuadPtsPerFace()*3;
  std::vector<Real> flux_vals_tmp(npts, 0), flux_vals_deriv_tmp(npts, 0);

  for (int i=0; i < npts; ++i)
  {
    flux_vals[i] = 0;
    flux_vals_deriv[i] = 0;
  }

  for (auto& bc : m_bcs)
  {
    bc->getValuedTair(face, t, sol_vals, flux_vals_tmp.data(), flux_vals_deriv_tmp.data());

    updateAndZero(flux_vals, flux_vals_tmp.data(), npts);
    updateAndZero(flux_vals_deriv, flux_vals_deriv_tmp.data(), npts);
  }
}

void CombinedAirWindSkyNeumannBC::getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar)
{
  int npts = getSurfDisc()->getNumQuadPtsPerFace();
  std::vector<Real> sol_vals_bar_tmp(npts, 0);

  for (int i=0; i < npts; ++i)
    sol_vals_bar[i] = 0;

  for (auto& bc : m_bcs)
  {
    bc->getValue_rev(face, t, sol_vals, sol_vals_bar_tmp.data(), flux_vals_bar);
    updateAndZero(sol_vals_bar, sol_vals_bar_tmp.data(), npts);
  }
}

void CombinedAirWindSkyNeumannBC::updateAndZero(Real* arr, Real* arr_tmp, int npts)
{
  for (int i=0; i < npts; ++i)
  {
#ifndef NDEBUG
    if (std::isnan(arr_tmp[i]))
      throw std::runtime_error("found nan");
#endif

    arr[i] += arr_tmp[i];
    arr_tmp[i] = 0;
  }
}


}  // namespace