#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#include "lagrange_basis.h"
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

    void computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs) override;
        
    void computeJacobian(DiscVectorPtr u, const Real t, linear_system::AssemblerPtr assembler) override;

    void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out) override;

    void computeMassMatrix(linear_system::AssemblerPtr assembler) override;

    void addVolumeGroupParams(const VolumeGroupParams& params) { m_params.push_back(params); }

    const VolumeGroupParams& getVolumeGroupParams(int idx) const { return m_params.at(idx); }

  private:
    void checkInitialization() override;  // check everything was initialized

    std::vector<VolumeGroupParams> m_params;
};


void applyDirichletValues(const HeatEquation& physics, const Real t, DiscVectorPtr u);

void applyMassMatrix(const HeatEquation& physics, DiscVectorPtr vec_in, DiscVectorPtr vec_out);

void applyMassMatrix(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const DofNumberingPtr dof_numbering,
                     const ArrayType<Real, 2>& arr_in, ArrayType<Real, 2>& arr_out);

void computeMassMatrix(const HeatEquation& physics, linear_system::AssemblerPtr assembler);

void computeMassMatrix(const VolDiscPtr vol_disc, const VolumeGroupParams& params, linear_system::AssemblerPtr assembler);

void computeVolumeTerm(const HeatEquation& physics, DiscVectorPtr u, DiscVectorPtr rhs);

void computeVolumeTerm(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                       const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);

void computeVolumeTerm2(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                        const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);

void computeVolumeJacobian(const HeatEquation& physics, DiscVectorPtr u, linear_system::AssemblerPtr assembler);

void computeVolumeTerm2Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                            linear_system::AssemblerPtr assembler);

void computeVolumeTerm3Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                           linear_system::AssemblerPtr assembler);

void computeSourceTerm(const HeatEquation& physics, Real t, DiscVectorPtr rhs);

void computeSourceTerm(const VolDiscPtr vol_disc, SourceTermPtr src, Real t,
                       ArrayType<Real, 2>& rhs_arr);

void computeNeumannBC(const HeatEquation& physics, const Real t, DiscVectorPtr u, DiscVectorPtr rhs);

void computeNeumannBC(NeumannBCPtr bc, DiscVectorPtr u, const Real t, DiscVectorPtr rhs);

void computeUnsteadyDirichletBC(const HeatEquation& physics, const Real t, DiscVectorPtr rhs);

void computeUnsteadyDirichletBC(DirichletBCPtr bc, const std::vector<VolDiscPtr>& vol_discs,
                                const std::vector<const VolumeGroupParams*>& vol_group_params,
                                const Real t, DiscVectorPtr rhs);

void printArray(DiscVectorPtr vec);

void printVector(DiscVectorPtr vec);



} // namespace

#endif