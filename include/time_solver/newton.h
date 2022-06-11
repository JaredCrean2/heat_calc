#ifndef TIME_SOLVER_NEWTON_H
#define TIME_SOLVER_NEWTON_H

#include <memory>
#include "ProjectDefs.h"
#include "time_solver/newton_result.h"

class DiscVector;
using DiscVectorPtr = std::shared_ptr<DiscVector>;

namespace linear_system {
  class LargeMatrix;
  using LargeMatrixPtr = std::shared_ptr<LargeMatrix>;
}

namespace timesolvers {

// solve a function f(u) = 0 using Newton's method
// The iteration is u_n+1 = u_n - (df/du_n)^-1 f(u_n)
class NewtonFunction
{
  public:
    virtual ~NewtonFunction() {}
    
    // function will be called before starting a new Newton solve.
    // This function facilitates using the same NewtonFunction object for
    // repeated solves
    virtual void resetForNewSolve() = 0;
  
    // compute f(u), overwriting f.  Computes norm of f if required
    virtual Real computeFunc(const DiscVectorPtr u, bool compute_norm, DiscVectorPtr f) = 0;

    // compute jac = df/du, overwriting jac
    virtual void computeJacobian(const DiscVectorPtr u, linear_system::LargeMatrixPtr jac) = 0;

    // create an empty vector
    virtual DiscVectorPtr createVector() = 0;
};

using NewtonFunctionPtr = std::shared_ptr<NewtonFunction>;


class NewtonSolver
{
  public:
    explicit NewtonSolver(NewtonFunctionPtr func, linear_system::LargeMatrixPtr jac);

    NewtonResult solve(DiscVectorPtr u, Real abs_tol, Real rel_tol, int itermax);

  private:
    void setupForSolve(DiscVectorPtr u, Real abs_tol, Real rel_tol, int itermax);

    void solveStep(DiscVectorPtr u);

    NewtonFunctionPtr m_func;
    linear_system::LargeMatrixPtr m_jac;
    DiscVectorPtr m_f;
    DiscVectorPtr m_delta_u;
    Real m_abs_tol;
    Real m_rel_tol;
    int m_itermax;
};

}

#endif