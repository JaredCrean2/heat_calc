#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#include "ProjectDefs.h"
#include "lagrange_basis.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/PhysicsModel.h"
#include "physics/heat/air_leakage.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/solar_position_calculator.h"
#include "quadrature.h"
#include <memory>

namespace Heat {

class VolumeGroupParams
{
  public:
    VolumeGroupParams (Real kappa=0, Real rho=0, Real Cp=0) :
      kappa(kappa),
      rho(rho),
      Cp(Cp)
    {}

    Real kappa; // thermal conductivity, W/(m K)
    Real rho;   // density (kg/m^2)
    Real Cp;    // specific heat capacityh (J/(kg K))
};

class InteriorAirTemperatureUpdator;

// solves heat equation rho Cp dT/dt = d/dx (kappa dT/dx) + S
// and S is a source term.  This module implements the
// two terms on the right hand side
// In weak form:
//  \int w rho Cp dT/dt = \int w d/dx (kappa dT/dx) + \int w S dOmega
//  \int w rho Cp dT/dt dOmega = - \int dw/dx kappa dT/dx dOmega + \int w kappa dT/dx n_i dGamma + \int w S dOmega = 0
//  \int w rho Cp dT/dt dOmega = - \int dw/dx kappa dT/dx dOmega + \int w h dGamma_h + \int w S dOmega = 0
// where h is the prescibes heat flux, Gamma is the boundary, and Gamma_h is the portion of the boundary
// where the heat flux is prescribed
// Note that w = 0 on Gamma_d, the portion of the boundary where the solution is specified, and
// Gamma_d \union Gamma_h = Gamma
class HeatEquation : public PhysicsModel
{
  public:
    HeatEquation(DiscPtr disc)
    : PhysicsModel(disc)
    {}

    void computeRhs(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, const Real t, DiscVectorPtr rhs) override;
        
    void computeJacobian(DiscVectorPtr u,  AuxiliaryEquationsStoragePtr u_aux, const Real t, linear_system::AssemblerPtr assembler) override;

    void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out) override;

    void computeMassMatrix(linear_system::AssemblerPtr assembler) override;

    void addVolumeGroupParams(const VolumeGroupParams& params) { m_params.push_back(params); }

    const VolumeGroupParams& getVolumeGroupParams(int idx) const { return m_params.at(idx); }

  private:
    void checkInitialization() override;  // check everything was initialized


    std::vector<VolumeGroupParams> m_params;
};


} // namespace

#endif
