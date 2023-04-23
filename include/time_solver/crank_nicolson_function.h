#ifndef TIME_SOLVER_CRANK_NICOLSON_FUNCTION
#define TIME_SOLVER_CRANK_NICOLSON_FUNCTION

#include "linear_system/augmented_assembler.h"
#include "physics/AuxiliaryEquations.h"
#include "time_solver/crank_nicolson_aux_equations.h"
#include "time_solver/newton.h"
#include "physics/PhysicsModel.h"

namespace timesolvers {

class CrankNicolsonAuxiliaryEquations;

class CrankNicolsonFunction : public NewtonFunction
{
  public:
    CrankNicolsonFunction(std::shared_ptr<PhysicsModel> physics_model, linear_system::LargeMatrixPtr mat,
                          Real t0, bool aux_eqns_combined_system=false);

    void resetForNewSolve() override ;
  
    Real computeFunc(const ArrayType<Real, 1>& u_np1, AuxiliaryEquationsStoragePtr u_aux_np1, bool compute_norm, ArrayType<Real, 1>& f_np1) override;

    // compute jac = df/du, overwriting jac
    void computeJacobian(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr jac) override;

    //void updateDependentQuantities(DiscVectorPtr u) override;

    virtual NewtonAuxiliaryEquationsPtr getAuxiliaryEquations() override { return m_aux_eqns; }

    void setTnp1(const ArrayType<Real, 1>& u_n, AuxiliaryEquationsStoragePtr u_aux_n, Real t_np1);

    // return the norm of the physics rhs (not the Crank Nicolson rhs) from
    // the last iteration when it was computed
    Real getLastPhysicsRhsNorm() { return m_last_physics_rhs_norm; }

  private:
    //void splitSolutionVector(const ArrayType<Real, 1>& combined_vec, ArrayType<Real, 1>& sol_vec,
    //                         AuxiliaryEquationsStoragePtr sol_aux);
//
    //void combineResidualVector(const ArrayType<Real, 1>& res_vec, AuxiliaryEquationsStoragePtr res_aux,
    //                            ArrayType<Real, 1>& res_combined);
    
    std::shared_ptr<PhysicsModel> m_physics_model;
    std::shared_ptr<CrankNicolsonAuxiliaryEquations> m_aux_eqns;
    linear_system::AssemblerPtr m_assembler;
    linear_system::AugmentedAssemblerPtr m_augmented_assembler;

    Real m_tn;
    bool m_aux_eqns_combined_system;
    //TODO: how many of these need to be true DiscVectors? Some of them may only
    //      require the vector part
    DiscVectorPtr m_un;
    AuxiliaryEquationsStoragePtr m_u_aux_n;
    DiscVectorPtr m_fn;
    AuxiliaryEquationsStoragePtr m_fn_aux;
    //AuxiliaryEquationsStoragePtr m_aux_tmp;

    Real m_tnp1;
    DiscVectorPtr m_delta_u;
    AuxiliaryEquationsStoragePtr m_delta_u_aux;
    DiscVectorPtr m_Mdelta_u;
    AuxiliaryEquationsStoragePtr m_Mdelta_u_aux;

    std::vector<DofInt> m_owned_dof_to_local;
    Real m_last_physics_rhs_norm = 0;
};

}

#endif