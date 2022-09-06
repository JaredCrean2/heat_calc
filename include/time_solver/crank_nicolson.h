#ifndef TIME_SOLVER_CRANK_NICOLSON
#define TIME_SOLVER_CRANK_NICOLSON

#include "ProjectDefs.h"
#include "discretization/disc_vector.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/PhysicsModel.h"
#include "linear_system/large_matrix.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix_factory.h"
#include "time_solver/newton.h"
#include "mesh/mesh.h"
#include "crank_nicolson_function.h"
#include "crank_nicolson_aux_equations.h"


namespace timesolvers {



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