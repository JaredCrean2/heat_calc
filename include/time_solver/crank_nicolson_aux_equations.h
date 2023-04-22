#ifndef TIME_SOLVER_CRANK_NICOLSON_AUX_EQNS
#define TIME_SOLVER_CRANK_NICOLSON_AUX_EQNS

#include "physics/AuxiliaryEquations.h"
#include "time_solver/newton.h"
#include "physics/PhysicsModel.h"


namespace timesolvers {

class CrankNicolsonAuxiliaryEquations : public NewtonAuxiliaryEquations
{
  public:
    CrankNicolsonAuxiliaryEquations(std::shared_ptr<PhysicsModel> physics_model, Real t0) :
      m_physics_model(physics_model),
      m_aux_eqns(physics_model->getAuxEquations()),
      m_tn(-1),
      m_tnp1(t0),
      m_un(makeDiscVector(physics_model->getDiscretization())),
      m_aux_un(makeAuxiliaryEquationsStorage(m_aux_eqns))
    {}

    virtual int getNumBlocks() const override { return m_aux_eqns->getNumBlocks(); }

    // returns the number of variables in the given block
    virtual int getBlockSize(int block) const override { return m_aux_eqns->getBlockSize(block); }

    virtual Real computeRhs(int block, const ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, bool compute_norm, ArrayType<Real, 1>& rhs) override;

    virtual void computeJacobian(int block, const ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr mat) override;

    virtual void multiplyOffDiagonal(int iblock, int jblock, const ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override;

    virtual AuxiliaryEquationsJacobiansPtr getJacobians() override
    {
      return m_aux_eqns->getJacobians();
    }

    virtual AuxiliaryEquationsStoragePtr createStorage() override;

  private:

    void setTnp1(DiscVectorPtr u_n, AuxiliaryEquationsStoragePtr u_aux_vec, Real t_np1);

    std::shared_ptr<PhysicsModel> m_physics_model;
    AuxiliaryEquationsPtr m_aux_eqns;
    Real m_tn;
    Real m_tnp1;
    DiscVectorPtr m_un;
    AuxiliaryEquationsStoragePtr m_aux_un;

    friend class CrankNicolsonFunction;
};

}

#endif