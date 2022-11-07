#ifndef TIME_SOLVER_CRANK_NICOLSON_FUNCTION
#define TIME_SOLVER_CRANK_NICOLSON_FUNCTION

#include "time_solver/crank_nicolson_aux_equations.h"
#include "time_solver/newton.h"
#include "physics/PhysicsModel.h"

namespace timesolvers {

class CrankNicolsonAuxiliaryEquations;

class CrankNicolsonFunction : public NewtonFunction
{
  public:
    CrankNicolsonFunction(std::shared_ptr<PhysicsModel> physics_model, linear_system::LargeMatrixPtr mat, Real t0);

    void resetForNewSolve() override ;
  
    Real computeFunc(const DiscVectorPtr u_np1, bool compute_norm, DiscVectorPtr f_np1) override;

    // compute jac = df/du, overwriting jac
    void computeJacobian(const DiscVectorPtr u, linear_system::LargeMatrixPtr jac) override;

    //void updateDependentQuantities(DiscVectorPtr u) override;

    virtual NewtonAuxiliaryEquationsPtr getAuxiliaryEquations() override { return m_aux_eqns; }

    // create an empty vector
    DiscVectorPtr createVector() override;

    void setTnp1(DiscVectorPtr u_n, Real t_np1);

  private:
    std::shared_ptr<PhysicsModel> m_physics_model;
    std::shared_ptr<CrankNicolsonAuxiliaryEquations> m_aux_eqns;
    linear_system::AssemblerPtr m_assembler;
    Real m_tn;
    //TODO: how many of these need to be true DiscVectors? Some of them may only
    //      require the vector part
    DiscVectorPtr m_un;
    DiscVectorPtr m_fn;

    Real m_tnp1;
    DiscVectorPtr m_delta_u;
    DiscVectorPtr m_Mdelta_u;
    std::vector<DofInt> m_owned_dof_to_local;
};

}

#endif