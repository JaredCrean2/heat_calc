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



}  // namespace