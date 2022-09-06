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



}  // namespace