#include "gtest/gtest.h"
#include "discretization/disc_vector.h"
#include "discretization/discretization.h"
#include "linear_system/assembler.h"
#include "physics/AuxiliaryEquations.h"
#include "time_solver/newton.h"
#include "mesh_helper.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern_mesh.h"
#include "time_solver/newton_result.h"

namespace {


linear_system::LargeMatrixOptsPetsc get_options()
{
  linear_system::LargeMatrixOptsPetsc opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;
  opts.petsc_opts["ksp_atol"] = "1e-14";
  opts.petsc_opts["ksp_rtol"] = "1e-50";
  //opts.petsc_opts["ksp_monitor"] = "";

  return opts;
}


class NewtonTester : public ::testing::Test,
                     public StandardDiscSetup
{
  protected:
    NewtonTester()
    {
      setup(5, 1, {false, false, false, false, false, false});

      auto opts = get_options();
      auto sparsity_pattern = std::make_shared<linear_system::SparsityPatternMesh>(mesh);
      mat = std::make_shared<linear_system::LargeMatrixPetsc>(mesh->getNumOwnedDofs(), mesh->getNumOwnedDofs(), opts, sparsity_pattern);
      u = makeDiscVector(disc);
    }


    std::shared_ptr<linear_system::LargeMatrixPetsc> mat;
    DiscVectorPtr u;
};

class NewtonTestAuxiliaryEquations : public timesolvers::NewtonAuxiliaryEquations
{
  public:
    NewtonTestAuxiliaryEquations(DiscPtr disc) :
      m_aux_eqns(std::make_shared<AuxiliaryEquationsNone>(disc))
    {}

    NewtonTestAuxiliaryEquations(DiscPtr disc, AuxiliaryEquationsPtr aux_eqns) :
      m_aux_eqns(aux_eqns)
    {}

    int getNumBlocks() const override { return m_aux_eqns->getNumBlocks(); }

    // returns the number of variables in the given block
    int getBlockSize(int block) const override { return m_aux_eqns->getBlockSize(block); }

    Real computeRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, bool compute_norm, ArrayType<Real, 1>& rhs) override 
    { 
      m_aux_eqns->computeRhs(block, u_vec, u_aux_vec, 0.0, rhs);

      if (compute_norm)
        return std::abs(rhs[0]);
      else
        return 0;
    }

    void computeJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr mat) override
    {
      auto assembler = std::make_shared<linear_system::SimpleAssembler>(mat);
      return m_aux_eqns->computeJacobian(block, u_vec, u_aux_vec, 0.0, assembler);
    }

    void multiplyOffDiagonal(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      m_aux_eqns->multiplyOffDiagonal(iblock, jblock, u_vec, u_aux_vec, 0.0, x, b);
    }

    AuxiliaryEquationsJacobiansPtr getJacobians() override
    {
      return m_aux_eqns->getJacobians();
    }

    AuxiliaryEquationsStoragePtr createStorage() override
    {
      return makeAuxiliaryEquationsStorage(m_aux_eqns);
    }

  private:
    AuxiliaryEquationsPtr m_aux_eqns;
};

// defines auxiliary system:
// [ a * x + c * y  = [0
// [ c * x + b * y] =  0]
class AuxiliaryEquationsCoupledBlock : public AuxiliaryEquations
{
  public:
    explicit AuxiliaryEquationsCoupledBlock(DiscPtr disc) :
      AuxiliaryEquations(disc)
    {
      m_jacs = std::make_shared<AuxiliaryEquationsJacobiansDense>(*this);
      //m_solutions = std::make_shared<AuxiliaryEquationStorage>(*this);
    }

    AuxiliaryEquationsJacobiansPtr getJacobians() override { return m_jacs; }

  protected:
    // return number of auxiliary sets of equations
    int getNumAuxiliaryBlocks() const override { return 2; }
    
    // returns number of variables in each block
    int getAuxiliaryBlockSize(int block) const override { return 1; }

    // each auxiliary block must be of the form du/dt = rhs(u, t).  This function computes the rhs
    void computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, ArrayType<Real, 1>& rhs) override
    {
      Real x = u_aux_vec->getVector(1)[0];
      Real y = u_aux_vec->getVector(2)[0];

      if (block == 0)
        rhs[0] = m_a * x + m_c * y - 3;
      else
        rhs[0] = m_c * x + m_b * y - 4;
    }

    void computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr mat) override
    {
      std::vector<DofInt> dofs = {0};
      ArrayType<Real, 2> jac(boost::extents[1][1]);
      jac[0][0] = 1;
      mat->assembleEntry(dofs, jac);
    }

    void multiplyAuxiliaryMassMatrix(int block, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      b[0] = x[0];
    }

    // compute the diagonal Jacobian block for the given block
    void computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, 
                                  Real t, linear_system::SimpleAssemblerPtr mat) override
    {
      std::vector<DofInt> dofs = {0};
      ArrayType<Real, 2> jac(boost::extents[1][1]);
      jac[0][0];

      if (block == 0)
        jac[0][0] = m_a;
      else
        jac[0][0] = m_b;  

      mat->assembleEntry(dofs, jac);
    }

    // compute the Jacobian-vector product for the block the couples the finite element problem to auxiliary block jblock
    void computeFiniteElementJacobianVectorProduct(int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec,
                                                   Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      for (int i=0; i < getBlockSize(0); ++i)
        b[i] = 0;
    }

    // compute the Jacobian-vector product for the block that couples auxiliary block i to auxiliary block j
    void computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec,
                                               Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      // the system we are solving is:
      // [a c  [x  = [b1
      //  c b]  y]    b2]
      // so the off-diagonal product is always c times the given value
      if (iblock == 0 || jblock == 0)
      {
        for (int i=0; i < b.shape()[0]; ++i)
          b[i] = 0;
      } else
        b[0] = m_c * x[0];
    }

  private:
    DiscPtr m_disc;
    AuxiliaryEquationsJacobiansPtr m_jacs;
    //AuxiliaryEquationsStoragePtr m_solutions;
    Real m_a = 1;
    Real m_b = 2;
    Real m_c = 0.1;

};


class NewtonTestFunc : public timesolvers::NewtonFunction
{
  public:
    explicit NewtonTestFunc(DiscPtr disc) : 
      m_disc(disc),
      m_aux_eqns(std::make_shared<NewtonTestAuxiliaryEquations>(disc))
    {
      disc->getMesh()->getOwnedLocalDofInfo(m_owned_dof_to_local);
      disc->getMesh()->getLocalToGlobalDofs(m_local_dof_to_global);
    }

    void setAuxiliaryEquations(timesolvers::NewtonAuxiliaryEquationsPtr aux_eqn)
    {
      m_aux_eqns = aux_eqn;
    }

    void resetForNewSolve() override {}
  
    Real computeFunc(const DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux_vec, bool compute_norm, DiscVectorPtr f) override
    {
      if (!u->isVectorCurrent())
        u->syncArrayToVector();

      f->set(0);
      auto& u_vec = u->getVector();
      auto& f_vec = f->getVector();

      for (int i=0; i < u->getNumDofs(); ++i)
        f_vec[i] = (u_vec[i] - 1)*(u_vec[i] - 1);

      double norm = 0.0, norm_global;
      if (compute_norm)
      {
        for (auto dof : m_owned_dof_to_local)
          norm += f_vec[dof]*f_vec[dof];

        MPI_Allreduce(&norm, &norm_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        norm = std::sqrt(norm_global);
      }

      return norm;
    }

    void computeJacobian(const DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr jac) override
    {    
      if (!u->isVectorCurrent())
        u->syncArrayToVector();

      std::vector<PetscInt> dofs(1);
      ArrayType<Real, 2> vals(boost::extents[1][1]);
      auto& u_vec = u->getVector();

      jac->zeroMatrix();
      for (int i=0; i < u->getNumDofs(); ++i)
      {
        dofs[0] = m_local_dof_to_global[i];
        vals[0][0] = 2*(u_vec[i] - 1);
        jac->assembleValues(dofs, vals);
      }
    }


    timesolvers::NewtonAuxiliaryEquationsPtr getAuxiliaryEquations() override { return m_aux_eqns; }

    virtual DiscVectorPtr createVector() override
    {
      return makeDiscVector(m_disc);
    }

  private:
    DiscPtr m_disc;
    std::vector<DofInt> m_owned_dof_to_local;
    std::vector<DofInt> m_local_dof_to_global;
    timesolvers::NewtonAuxiliaryEquationsPtr m_aux_eqns;
};


}  // namespace


TEST_F(NewtonTester, Quadratic)
{
  auto func = std::make_shared<NewtonTestFunc>(disc);
  const Real abstol = 1e-13;
  const int itermax = 100;

  timesolvers::NewtonSolver newton(func, mat);
  u->set(2);
  auto u_aux = func->getAuxiliaryEquations()->createStorage();
  timesolvers::NewtonOpts opts;
  opts.nonlinear_abs_tol = abstol;
  opts.nonlinear_itermax = itermax;
  opts.nonlinear_rel_tol = -1;
  opts.linear_itermax = 5;
  auto result = newton.solve(u, u_aux, opts);

  if (!u->isVectorCurrent())
    u->syncArrayToVector();

  auto& u_vec = u->getVector();
  for (int i=0; i < u->getNumDofs(); ++i)
    EXPECT_NEAR(u_vec[i], 1.0, std::sqrt(abstol));

  EXPECT_TRUE(result.isConverged());
  //EXPECT_EQ(result.getNIters(), 1);

  EXPECT_TRUE(result.isAbstolSatisfied());
  EXPECT_FALSE(result.isReltolSatisfied());
  EXPECT_TRUE(result.isItermaxSatisfied());

  EXPECT_EQ(result.getAbsTol(), abstol);
  EXPECT_EQ(result.getRelTol(), -1);
}

TEST_F(NewtonTester, CoupledAuxiliaryBlock)
{
  auto func            = std::make_shared<NewtonTestFunc>(disc);
  auto aux_eqns        = std::make_shared<AuxiliaryEquationsCoupledBlock>(disc);
  auto newton_aux_eqns = std::make_shared<NewtonTestAuxiliaryEquations>(disc, aux_eqns);
  auto u_aux           = makeAuxiliaryEquationsStorage(aux_eqns);
  func->setAuxiliaryEquations(newton_aux_eqns);

  const Real abstol = 1e-13;
  const int itermax = 100;

  timesolvers::NewtonSolver newton(func, mat);
  u->set(2);
  u_aux->getVector(1)[0] = 1;
  u_aux->getVector(2)[0] = 2;

  timesolvers::NewtonOpts opts;
  opts.nonlinear_abs_tol = abstol;
  opts.nonlinear_itermax = itermax;
  opts.nonlinear_rel_tol = -1;
  auto result = newton.solve(u, u_aux, opts);

  if (!u->isVectorCurrent())
    u->syncArrayToVector();

  auto& u_vec = u->getVector();
  for (int i=0; i < u->getNumDofs(); ++i)
    EXPECT_NEAR(u_vec[i], 1.0, std::sqrt(abstol));

  EXPECT_TRUE(result.isConverged());
  //EXPECT_EQ(result.getNIters(), 1);

  EXPECT_TRUE(result.isAbstolSatisfied());
  EXPECT_FALSE(result.isReltolSatisfied());
  EXPECT_TRUE(result.isItermaxSatisfied());

  EXPECT_EQ(result.getAbsTol(), abstol);
  EXPECT_EQ(result.getRelTol(), -1);

  Real x = u_aux->getVector(1)[0];
  Real y = u_aux->getVector(2)[0];

  EXPECT_NEAR(x, 2.81407035175875, abstol);
  EXPECT_NEAR(y, 1.859296482412063, abstol);
}