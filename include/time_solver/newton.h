#ifndef TIME_SOLVER_NEWTON_H
#define TIME_SOLVER_NEWTON_H

#include <memory>
#include "ProjectDefs.h"
#include "physics/AuxiliaryEquations.h"
#include "time_solver/newton_result.h"

class DiscVector;
using DiscVectorPtr = std::shared_ptr<DiscVector>;

namespace linear_system {
  class LargeMatrix;
  using LargeMatrixPtr = std::shared_ptr<LargeMatrix>;
}

namespace timesolvers {

class NewtonAuxiliaryEquations
{
  public:
    virtual ~NewtonAuxiliaryEquations() {};

    virtual int getNumBlocks() const = 0;

    // returns the number of variables in the given block
    virtual int getBlockSize(int block) const = 0;

    virtual Real computeRhs(int block, DiscVectorPtr u_vec, bool compute_norm, ArrayType<Real, 1>& rhs) = 0;

    virtual void computeJacobian(int block, DiscVectorPtr u_vec, linear_system::LargeMatrixPtr mat) =0;

    virtual void multiplyOffDiagonal(int iblock, int jblock, DiscVectorPtr u_vec, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) = 0;

    virtual void setBlockSolution(int block, const ArrayType<Real, 1>& vals) = 0;

    virtual ArrayType<Real, 1>& getBlockSolution(int block) = 0;

    virtual AuxiliaryEquationsJacobiansPtr getJacobians() = 0;

    virtual AuxiliaryEquationsStoragePtr createStorage() = 0;
};

using NewtonAuxiliaryEquationsPtr = std::shared_ptr<NewtonAuxiliaryEquations>;


// solve a function f(u) = 0 using Newton's method
// The iteration is u_n+1 = u_n - (df/du_n)^-1 f(u_n)
// If there are auxiliary equations, block Gauss-Seidel is used
// to converge them too.
class NewtonFunction
{
  public:
    virtual ~NewtonFunction() {}
    
    // function will be called before starting a new Newton solve.
    // This function facilitates using the same NewtonFunction object for
    // repeated solves
    virtual void resetForNewSolve() = 0;
  
    //TODO: change the API: the completeTimestep function should be the only
    //      function that changes the state u at which the function and
    //      jacobian are evaluated.  Remove u as an explicit argument of
    //      computeFunc, computeJacobian

    // compute f(u), overwriting f.  Computes norm of f if required
    virtual Real computeFunc(const DiscVectorPtr u, bool compute_norm, DiscVectorPtr f) = 0;

    // compute jac = df/du, overwriting jac
    virtual void computeJacobian(const DiscVectorPtr u, linear_system::LargeMatrixPtr jac) = 0;

    virtual NewtonAuxiliaryEquationsPtr getAuxiliaryEquations() = 0;

    //virtual void updateDependentQuantities(DiscVectorPtr u) {};

    // create an empty vector
    virtual DiscVectorPtr createVector() = 0;
};

using NewtonFunctionPtr = std::shared_ptr<NewtonFunction>;


class NewtonSolver
{
  public:
    explicit NewtonSolver(NewtonFunctionPtr func, linear_system::LargeMatrixPtr jac);

    NewtonResult solve(DiscVectorPtr u, NewtonOpts opts);

  private:
    void setupForSolve(DiscVectorPtr u, NewtonOpts opts);

    void solveStep(DiscVectorPtr u);

    Real computeRhsAndNorm(DiscVectorPtr u);

    void computeJacobians(DiscVectorPtr u);

    // returns some kind of relative norm of delta_u
    Real gaussSeidelStep(DiscVectorPtr u);

    void computeLinearResidual(DiscVectorPtr u_vec);


    NewtonFunctionPtr m_func;
    NewtonOpts m_opts;
    linear_system::LargeMatrixPtr m_jac;
    AuxiliaryEquationsJacobiansPtr m_aux_jacs;
    DiscVectorPtr m_f;
    DiscVectorPtr m_delta_u;
    AuxiliaryEquationsStoragePtr m_aux_u;
    AuxiliaryEquationsStoragePtr m_aux_delta_u;
    AuxiliaryEquationsStoragePtr m_aux_rhs;
};

}

#endif