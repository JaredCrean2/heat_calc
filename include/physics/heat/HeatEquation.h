#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#include "lagrange_basis.h"
#include "mesh/mesh.h"
#include "physics/PhysicsModel.h"
#include "quadrature.h"

namespace Heat {

class VolumeGroupParams
{
  public:
    VolumeGroupParams (Real kappa, Real rho, Real Cp) :
      kappa(kappa),
      rho(rho),
      Cp(Cp)
    {}

    Real kappa; // thermal conductivity, W/(m K)
    Real rho;   // density (kg/m^2)
    Real Cp;    // specific heat capacityh (J/(kg K))
};

// solves heat equation dT/dt = d/dx (alpha dT/dx) + S
// where alpha = kappa / (rho Cp)
// and S is a source term
// To allow nonlinear problems to be solves, the equation is implemented as
// dT/dt - d/dx (alpha dT/dx) - S = 0, where this module only evaluates
// the final two terms on the left hand side.
// In weak form:
//  \int w dT/dt - \int w d/dx (alpha dT/dx) - \int w S dOmega = 0
//  \int w dT/dt dOmega + \int dw/dx alpha dT/dx dOmega - \int w alpha dT/dx n_i dGamma + - \int w S dOmega = 0
//  \int w dT/dt dOmega + \int dw/dx alpha dT/dx dOmega - \int w h dGamma_h + - \int w S dOmega = 0
// where h is the prescibes heat flux, Gamma is the boundary, and Gamma_h is the portion of the boundary
// where the heat flux is prescribed
class HeatEquation : public PhysicsModel
{
  public:
    HeatEquation(DiscPtr disc)
    : PhysicsModel(disc)
    {}

    void computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs) override;

    void addVolumeGroupParams(const VolumeGroupParams& params) { m_params.push_back(params); }

    const VolumeGroupParams& getVolumeGroupParams(int idx) const { return m_params.at(idx); }

  private:
    void checkInitialization() override;  // check everything was initialized

    std::vector<VolumeGroupParams> m_params;
};


void applyDirichletValues(const HeatEquation& physics, const Real t, DiscVectorPtr u);

void computeVolumeTerm(const HeatEquation& physics, DiscVectorPtr u, DiscVectorPtr rhs);

void computeVolumeTerm(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                       const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);

void computeSourceTerm(const HeatEquation& physics, Real t, DiscVectorPtr rhs);

void computeSourceTerm(const VolDiscPtr vol_disc, SourceTermPtr src, Real t,
                       ArrayType<Real, 2>& rhs_arr);

} // namespace

#endif