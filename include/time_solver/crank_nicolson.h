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


namespace timesolvers {

class CrankNicolsonAuxiliaryEquations : public NewtonAuxiliaryEquations
{
  public:
    CrankNicolsonAuxiliaryEquations(std::shared_ptr<PhysicsModel> physics_model, Real t0) :
      m_aux_eqns(physics_model->getAuxEquations()),
      m_tn(-1),
      m_tnp1(t0),
      m_un(makeDiscVector(physics_model->getDiscretization())),
      m_aux_un(*m_aux_eqns),
      m_aux_unp1(*m_aux_eqns)
    {}

    virtual int getNumBlocks() const override { return m_aux_eqns->getNumBlocks(); }

    // returns the number of variables in the given block
    virtual int getBlockSize(int block) const override { return m_aux_eqns->getBlockSize(block); }

    virtual void computeRhs(int block, DiscVectorPtr u_vec, ArrayType<Real, 1>& rhs) override
    {
      int num_vars = getBlockSize(block);

      // compute M * (u_np1 - u_n)/delta_t
      ArrayType<Real, 1> delta_u(boost::extents[num_vars]);
      auto& u_np1 = m_aux_unp1.getVector(block);
      auto& u_n   = m_aux_un.getVector(block);
      for (int i=0; i < num_vars; ++i)
        delta_u[i] = (u_np1[i] - u_n[i])/(m_tnp1 - m_tn);
      m_aux_eqns->multiplyMassMatrix(block, m_tnp1, delta_u, rhs);

      // compute 1/2(f(u_np1, t_np1) + f(u_n, t_n))
      ArrayType<Real, 1> rhs_tmp(boost::extents[num_vars]), rhs_tmp2(boost::extents[num_vars]);
      m_aux_eqns->computeRhs(block, m_un, m_tn, rhs_tmp);
      m_aux_eqns->computeRhs(block, u_vec, m_tnp1, rhs_tmp2);

      for (int i=0; i < num_vars; ++i)
        rhs[i] -= 0.5*(rhs_tmp[i] + rhs_tmp2[i]);
    }

    virtual void computeJacobian(int block, DiscVectorPtr u_vec, linear_system::LargeMatrixPtr mat) override
    {
      auto assembler = std::make_shared<linear_system::SimpleAssembler>(mat);
      assembler->setAlpha(1.0/(m_tnp1 - m_tn));
      m_aux_eqns->computeMassMatrix(block, m_tnp1, assembler);

      assembler->setAlpha(0.5);
      m_aux_eqns->computeJacobian(block, u_vec, m_tnp1, assembler);
    }

    virtual void multiplyOffDiagonal(int iblock, int jblock, DiscVectorPtr u_vec, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      m_aux_eqns->multiplyOffDiagonal(iblock, jblock, u_vec, m_tnp1, x, b);
      int num_vars = getBlockSize(iblock);
      for (int i=0; i < num_vars; ++i)
        b[i] *= -0.5;
    }

    virtual void setBlockSolution(int block, const ArrayType<Real, 1>& vals) override
    {
      m_aux_eqns->setBlockSolution(block, vals);
    }

    virtual AuxiliaryEquationsJacobiansPtr getJacobians() override
    {
      return m_aux_eqns->getJacobians();
    }

    virtual AuxiliaryEquationsStoragePtr createStorage() override;


  private:

    void setTnp1(DiscVectorPtr u_n, Real t_np1)
    {
      *m_un  = *u_n;
      m_tn   = m_tnp1;
      m_tnp1 = t_np1;
    }

    AuxiliaryEquationsPtr m_aux_eqns;
    Real m_tn;
    Real m_tnp1;
    DiscVectorPtr m_un;
    AuxiliaryEquationStorage m_aux_un;
    AuxiliaryEquationStorage m_aux_unp1;

    friend class CrankNicolsonFunction;
};


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