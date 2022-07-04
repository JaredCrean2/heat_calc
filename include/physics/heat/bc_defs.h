#ifndef HEAT_BC_DEFS_H
#define HEAT_BC_DEFS_H

#include "discretization/NeumannBC.h"
#include "discretization/surface_discretization.h"

namespace Heat {

class NewtonCooling : public NeumannBC
{
  public:
    NewtonCooling(SurfDiscPtr surf, Real heat_transfer_coeff) :
      NeumannBC(surf, true),
      m_heat_transfer_coeff(heat_transfer_coeff)
    {}

    void setExternalTemperature(Real temp) { m_temp = temp; }

    Real getExternalTemperature() const { return m_temp; }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
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


    void getValueDeriv(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals_deriv)
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

  private:
    Real m_heat_transfer_coeff; // Units W/(m^2 K)
    Real m_temp = std::numeric_limits<Real>::min();
};


// TARP (Thermal Analysis Research Program) model.
// See Section 3.5.5.2 of EnergyPlus Engineering reference Version 9.5
//class Tarp : public NeumannBC
//{
//  public:
//    Tarp(SurfDiscPtr surf) :
//      NeumannBC(surf, true)
//    {}
//
//    void setAirTemperature(Real temp) { m_air_temp = temp; }
//
//    Real getAirTemperature() const { return m_air_temp; }
//
//    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals)
//    {
//      auto& normals = m_surf->normals;
//      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
//      {
//        Real nx = normals[face][i][0];
//        Real ny = normals[face][i][1];
//        Real nz = normals[face][i][2];
//        Real mag     = std::sqrt(nx*nx + ny*ny + nz*nz);
//        std::array<Real, 3> unit_normal{nx/mag, ny/mag, nz/mag};
//        Real delta_t = m_air_temp - sol_vals[i];
//        Real heat_transfer_coeff = computeHeatTransferCoeff(delta_t, unit_normal);
//        Real term    = heat_transfer_coeff * delta_t / mag;
//
//        flux_vals[i]                                      = nx * term;
//        flux_vals[i + m_surf->getNumQuadPtsPerFace()]     = ny * term;
//        flux_vals[i + 2 * m_surf->getNumQuadPtsPerFace()] = nz * term;
//      }
//    }
//
//
//    void getValueDeriv(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals_deriv)
//    {
//      auto& normals = m_surf->normals;
//      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
//      {
//        Real nx = normals[face][i][0];
//        Real ny = normals[face][i][1];
//        Real nz = normals[face][i][2];
//        Real mag     = std::sqrt(nx*nx + ny*ny + nz*nz);
//        std::array<Real, 3> unit_normal{nx/mag, ny/mag, nz/mag};
//        Real delta_t = m_air_temp - sol_vals[i];
//        Real heat_transfer_coeff_dot = 0.0;
//        /*Real heat_transfer_coeff =*/ computeHeatTransferCoeffDeriv(delta_t, unit_normal, heat_transfer_coeff_dot);
//        //Real term    = heat_transfer_coeff * delta_t / mag;
//        Real term_dot    = heat_transfer_coeff_dot * delta_t / mag;
//
//        flux_vals_deriv[i]                                      = nx * term_dot;
//        flux_vals_deriv[i + m_surf->getNumQuadPtsPerFace()]     = ny * term_dot;
//        flux_vals_deriv[i + 2 * m_surf->getNumQuadPtsPerFace()] = nz * term_dot;
//      }
//    }
//
//  private:
//
//
//
//
//};
//
}

#endif