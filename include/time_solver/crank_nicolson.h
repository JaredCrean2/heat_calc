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
#include "time_solver/timestep_controller.h"
#include "time_stepper_opts.h"


namespace timesolvers {

class CrankNicolson
{
  public:
    CrankNicolson(std::shared_ptr<PhysicsModel> physics_model, DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, TimeStepperOpts opts);

    void solve();

  private:
    void advanceTimestep(Real t_new, Real delta_t);

    int numWholeSteps();

    double finalStepSize();

    std::shared_ptr<PhysicsModel> m_physics_model;
    std::shared_ptr<AuxiliaryEquations> m_aux_eqns;
    linear_system::LargeMatrixPtr m_matrix;
    DiscVectorPtr m_u;
    AuxiliaryEquationsStoragePtr m_u_aux;
    TimeStepperOpts m_opts;
    std::shared_ptr<CrankNicolsonFunction> m_func;
    std::shared_ptr<NewtonSolver> m_newton;
};

} // namespace

#endif