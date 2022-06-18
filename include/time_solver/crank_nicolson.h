#ifndef TIME_SOLVER_CRANK_NICOLSON
#define TIME_SOLVER_CRANK_NICOLSON

#include "ProjectDefs.h"
#include "discretization/disc_vector.h"
#include "physics/PhysicsModel.h"
#include "linear_system/large_matrix.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix_factory.h"
#include "time_solver/newton.h"
#include "mesh/mesh.h"


namespace timesolvers {

class CrankNicolsonFunction : public NewtonFunction
{
  public:
    CrankNicolsonFunction(std::shared_ptr<PhysicsModel> physics_model, linear_system::LargeMatrixPtr mat, Real t0);

    void resetForNewSolve() override ;
  
    Real computeFunc(const DiscVectorPtr u_np1, bool compute_norm, DiscVectorPtr f_np1) override;

    // compute jac = df/du, overwriting jac
    void computeJacobian(const DiscVectorPtr u, linear_system::LargeMatrixPtr jac) override;

    // create an empty vector
    DiscVectorPtr createVector() override;

    void setTnp1(DiscVectorPtr u_n, Real t_np1);

  private:
    std::shared_ptr<PhysicsModel> m_physics_model;
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


struct TimeStepperOpts
{
  Real t_start = 0;
  Real t_end   = 0;
  Real delta_t = 0;
  linear_system::LargeMatrixType mat_type = linear_system::LargeMatrixType::Unknown;
  std::shared_ptr<linear_system::LargeMatrixOpts> matrix_opts = nullptr;
  Real nonlinear_abs_tol = -1;
  Real nonlinear_rel_tol = -1;
  int nonlinear_itermax  = -1;
};

void checkTimeStepperOpts(const TimeStepperOpts& opts, bool check_implicit=true);


class CrankNicolson
{
  public:
    CrankNicolson(std::shared_ptr<PhysicsModel> physics_model, DiscVectorPtr u, TimeStepperOpts opts);

    void solve();

  private:
    void advanceTimestep(Real t_new, Real delta_t);

    int numWholeSteps();

    double finalStepSize();

    std::shared_ptr<PhysicsModel> m_physics_model;
    linear_system::LargeMatrixPtr m_matrix;
    DiscVectorPtr m_u;
    TimeStepperOpts m_opts;
    std::shared_ptr<CrankNicolsonFunction> m_func;
    std::shared_ptr<NewtonSolver> m_newton;
};

} // namespace

#endif