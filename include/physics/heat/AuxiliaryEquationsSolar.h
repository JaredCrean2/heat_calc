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
    {
      m_solutions = std::make_shared<AuxiliaryEquationStorage>(*this);
    }

  protected:

    // return number of auxiliary sets of equations
    int getNumAuxiliaryBlocks() const override { return 1; }
    
    // returns number of variables in each block
    int getAuxiliaryBlockSize(int block) const override { return 1; }

    // each auxiliary block must be of the form M du/dt = rhs(u, t).  This function computes the rhs
    void computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, Real t, ArrayType<Real, 1>& rhs) override
    {
      Real interior_temp = getAuxiliaryBlockSolution(block)[0];
      rhs[0] = m_air_temp->computeNetFlux(u_vec, interior_temp, t);
    }

    void computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr assembler) override
    {
      std::vector<DofInt> dofs = {0};
      ArrayType<Real, 2> mat(boost::extents[1][1]);
      mat[0][0] = m_air_temp->getThermalMass();
      assembler->assembleEntry(dofs, mat);
    }

    void multiplyAuxiliaryMassMatrix(int block, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      Real val = m_air_temp->getThermalMass();
      b[0] = x[0] * val;
    }

    // compute the diagonal Jacobian block for the given block
    void computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, Real t, linear_system::SimpleAssemblerPtr mat) override
    {
      Real interior_temp = getAuxiliaryBlockSolution(0)[0];
      Real val = m_air_temp->computeNetFluxJacobian(u_vec, interior_temp, t);

      std::vector<DofInt> dofs = {0};
      ArrayType<Real, 2> vals(boost::extents[1][1]);
      vals[0][0] = val;
      mat->assembleEntry(dofs, vals);
    }

    // compute the diagonal Jacobian block for the given block
    //virtual void computeAuxiliaryJacobian(int block, Real t, linear_system::SimpleAssemblerPtr mat) = 0;

    // compute the Jacobian-vector product for the block the couples the finite element problem to auxiliary block jblock
    void computeFiniteElementJacobianVectorProduct(int jblock, DiscVectorPtr u_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      m_heat_eqn.computedRdTinterior_airProduct(u_vec, t, x[0], b);
    }

    // compute the Jacobian-vector product for the block that couples auxiliary block i to auxiliary block j
    void computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, DiscVectorPtr u_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      Real t_interior = getAuxiliaryBlockSolution(0)[0];
      auto u_bar = makeDiscVector(m_heat_eqn.getDiscretization());  //TODO: cache this
      u_bar->set(0);
      m_air_temp->computeNetFlux_rev(u_vec, t_interior, t, u_bar, 1);

      if (!u_bar->isVectorCurrent())
        u_bar->syncArrayToVector();

      auto& u_bar_vec = u_bar->getVector();
      b[0] = 0;
      for (size_t i=0; i < u_bar_vec.shape()[0]; ++i)
        b[0] += u_bar_vec[i] * x[i];
    }

    AuxiliaryEquationsJacobiansPtr getJacobians() override { return m_jacs; }

    void setAuxiliaryBlockSolution(int block, const ArrayType<Real, 1>& vals) override;

    ArrayType<Real, 1>& getAuxiliaryBlockSolution(int block) override;


  private:
    HeatEquationSolar& m_heat_eqn;
    std::shared_ptr<InteriorAirTemperatureUpdator> m_air_temp;
    AuxiliaryEquationsJacobiansPtr m_jacs;
    AuxiliaryEquationsStoragePtr m_solutions;
};

}

#endif