#ifndef PHYSICS_HEAT_AUXILIARY_EQNS_SOLAR_H
#define PHYSICS_HEAT_AUXILIARY_EQNS_SOLAR_H

#include "HeatEquationSolar.h"
#include "interior_temperature_update.h"
#include "physics/AuxiliaryEquations.h"

namespace Heat {

class AuxiliaryEquationsSolar : public AuxiliaryEquations
{
  public:
    explicit AuxiliaryEquationsSolar(HeatEquationSolar& heat_eqn, std::shared_ptr<InteriorAirTemperatureUpdator> air_temp) :
      AuxiliaryEquations(heat_eqn.getDiscretization()),
      m_heat_eqn(heat_eqn),
      m_air_temp(air_temp),
      m_jacs(std::make_shared<AuxiliaryEquationsJacobiansDense>(*this))
    {}

  protected:

    // return number of auxiliary sets of equations
    int getNumAuxiliaryBlocks() const override { return 1; }
    
    // returns number of variables in each block
    int getAuxiliaryBlockSize(int block) const override { return 1; }

    // each auxiliary block must be of the form M du/dt = rhs(u, t).  This function computes the rhs
    void computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, ArrayType<Real, 1>& rhs) override;

    void computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr assembler) override;

    void multiplyAuxiliaryMassMatrix(int block, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override;

    // compute the diagonal Jacobian block for the given block
    void computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, linear_system::SimpleAssemblerPtr mat) override;

    // compute the diagonal Jacobian block for the given block
    //virtual void computeAuxiliaryJacobian(int block, Real t, linear_system::SimpleAssemblerPtr mat) = 0;

    // compute the Jacobian-vector product for the block the couples the finite element problem to auxiliary block jblock
    void computeFiniteElementJacobianVectorProduct(int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, 
                                                   const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override;

    // compute the Jacobian-vector product for the block that couples auxiliary block i to auxiliary block j
    void computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, 
                                                const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override;

    AuxiliaryEquationsJacobiansPtr getJacobians() override { return m_jacs; }

    //-------------------------------------------------------------------------
    // For augmented system
    void computeAuxiliaryJacobianDiagonalBlock(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec,
                                               Real t, linear_system::AugmentedAssemblerPtr mat) override;

    // compute the block that couples the finite element jacobian to the jth auxiliary block
    void computeFiniteElementJacobianOffDiagonallBlock(int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec,
                                                       Real t, linear_system::AugmentedAssemblerPtr mat) override;
    // assembles block that couples iblock to jblock
    void computeAuxiliaryJacobianOffDiagonalBlock(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec,
                                                  Real t, linear_system::AugmentedAssemblerPtr mat) override;

  private:
    HeatEquationSolar& m_heat_eqn;
    std::shared_ptr<InteriorAirTemperatureUpdator> m_air_temp;
    AuxiliaryEquationsJacobiansPtr m_jacs;
};

}

#endif