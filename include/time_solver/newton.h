#ifndef TIME_SOLVER_NEWTON_H
#define TIME_SOLVER_NEWTON_H

#include <memory>
#include "ProjectDefs.h"
#include "physics/AuxiliaryEquations.h"
#include "time_solver/newton_result.h"


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

    virtual Real computeRhs(int block, const ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, 
                            bool compute_norm, ArrayType<Real, 1>& rhs) = 0;

    virtual void computeJacobian(int block, const ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec,
                                 linear_system::LargeMatrixPtr mat) =0;

    virtual void multiplyOffDiagonal(int iblock, int jblock, const ArrayType<Real, 1>& u_vec, 
                                     AuxiliaryEquationsStoragePtr u_aux_vec, 
                                     const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) = 0;

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


    // compute f(u), overwriting f.  Computes norm of f if required
    virtual Real computeFunc(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, bool compute_norm, ArrayType<Real, 1>& f) = 0;

    // compute jac = df/du, overwriting jac
    virtual void computeJacobian(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr jac) = 0;

    virtual NewtonAuxiliaryEquationsPtr getAuxiliaryEquations() = 0;
};

using NewtonFunctionPtr = std::shared_ptr<NewtonFunction>;


class NewtonSolver
{
  public:
    explicit NewtonSolver(NewtonFunctionPtr func, linear_system::LargeMatrixPtr jac, MPI_Comm comm=MPI_COMM_WORLD);

    NewtonResult solve(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, NewtonOpts opts);

  private:
    void setupForSolve(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, NewtonOpts opts);

    void solveStep(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    void updateNonlinearSolution(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    Real computeRhsAndNorm(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    void computeJacobians(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    void checkJacobianFiniteDifference(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    // returns some kind of relative norm of delta_u
    Real gaussSeidelStep(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    Real gaussSeidelStepFirstRow(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    Real gaussSeidelStepOtherRows(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec);

    void gaussSeidelComputeRhs(int iblock, const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, ArrayType<Real, 1>& rhs);

    Real updateLinearSolution(const ArrayType<Real, 1>& delta_u_tmp, ArrayType<Real, 1>& delta_u);

    void computeLinearResidual(ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec);


    NewtonFunctionPtr m_func;
    NewtonOpts m_opts;
    linear_system::LargeMatrixPtr m_jac;
    AuxiliaryEquationsJacobiansPtr m_aux_jacs;
    ArrayType<Real, 1> m_f;
    ArrayType<Real, 1> m_delta_u;

    // linear solver temporaries
    AuxiliaryEquationsStoragePtr m_aux_delta_u;
    AuxiliaryEquationsStoragePtr m_aux_rhs;
    //MPI_Comm m_comm;
    bool m_am_i_root;
};

}

#endif