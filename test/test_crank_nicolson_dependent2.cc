#include "gtest/gtest.h"
#include "mesh_helper.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/PhysicsModel.h"
#include "time_solver/crank_nicolson.h"

namespace {

// one block is a function only of t, and the second block is a
// function of t and the other block.  This enum determines
// which block dependends on the other
enum class BlockDependent
{
  FirstBlock,
  SecondBlock
};

class CNDependentTester2 : public StandardDiscSetup,
                           public testing::Test
{
  protected:
    CNDependentTester2()
    {
      Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 3, 3, 3);
      setup(3, 1, spec, {false, false, false, false, false, false});
      //setup(3, 1, spec, {true, true, true, true, true, true});
    }
};


class AuxiliaryEquationsTest : public AuxiliaryEquations
{
  public:
    explicit AuxiliaryEquationsTest(DiscPtr disc, BlockDependent dependent) :
      AuxiliaryEquations(disc),
      m_disc(disc),
      m_jacs(std::make_shared<AuxiliaryEquationsJacobiansDense>(*this)),
      m_block_dependent(dependent)
    {}


    virtual AuxiliaryEquationsJacobiansPtr getJacobians() override { return m_jacs; }

  protected:
    // return number of auxiliary sets of equations
    int getNumAuxiliaryBlocks() const override { return 1; }
    
    // returns number of variables in each block
    int getAuxiliaryBlockSize(int block) const override { return m_disc->getDofNumbering()->getNumOwnedDofs(); }

    // u_aux(t) = t^2 + 0.1*t + 1, d u_aux/dt = 2*t + 0.1*t = 2*t + 0.2*sqrt(u)
    void computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, ArrayType<Real, 1>& rhs) override
    {
      auto& u_vec_vec = u_vec->getVector();
      for (int i=0; i < rhs.shape()[0]; ++i)
        if (m_block_dependent == BlockDependent::FirstBlock)
        {
          rhs[i] = 2*t;
        } else
        {
          rhs[i] = 2*t + 0.2*std::sqrt(u_vec_vec[i]);
        }
    }

    void computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr mat) override
    {
      std::vector<DofInt> dofs(1);
      ArrayType<Real, 2> jac(boost::extents[1][1]);
      for (int i=0; i < getAuxiliaryBlockSize(block); ++i)
      {
        dofs[0] = i;
        jac[0][0] = 1;
        mat->assembleEntry(dofs, jac);
      }
    }

    void multiplyAuxiliaryMassMatrix(int block, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      for (int i=0; i < x.shape()[0]; ++i)
        b[i] = x[i];
    }

    // compute the diagonal Jacobian block for the given block
    void computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, linear_system::SimpleAssemblerPtr mat) override
    {
      // nothing to do here: the diagonal block of the Jacobian is 0
    }

    // compute the diagonal Jacobian block for the given block
    //virtual void computeAuxiliaryJacobian(int block, Real t, linear_system::SimpleAssemblerPtr mat) = 0;

    // compute the Jacobian-vector product for the block the couples the finite element problem to auxiliary block jblock
    void computeFiniteElementJacobianVectorProduct(int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t,
                                                   const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      auto& u_aux_vec_vec = u_aux_vec->getVector(1);
      for (int i=0; i < b.shape()[0]; ++i)
        if (m_block_dependent == BlockDependent::FirstBlock)
        {
          b[i] = 0.1/std::sqrt(u_aux_vec_vec[i]) * x[i];
        } else
        {
          b[i] = 0;
        }
    }

    // compute the Jacobian-vector product for the block that couples auxiliary block i to auxiliary block j
    void computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t,
                                               const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      if (iblock == 1 && jblock == 0 && m_block_dependent == BlockDependent::SecondBlock)
      {
        auto& u_vec_vec = u_vec->getVector();
        for (int i=0; i < b.shape()[0]; ++i)
          b[i] = 0.1/std::sqrt(u_vec_vec[i]) * x[i];
      } else
      {
        for (int i=0; i < b.shape()[0]; ++i)
          b[i] = 0;
      }
    }

  private:
    DiscPtr m_disc;
    AuxiliaryEquationsJacobiansPtr m_jacs;
    BlockDependent m_block_dependent;
};

class PhysicsModelTest : public PhysicsModel
{
  public:
    explicit PhysicsModelTest(DiscPtr disc, BlockDependent dependent) :
        PhysicsModel(disc),
        m_aux_eqns(std::make_shared<AuxiliaryEquationsTest>(disc, dependent)),
        m_block_dependent(dependent)
    {}

    // u(t) = t^2, therefore du/dt = 2*t
    void computeRhs(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, const Real t, DiscVectorPtr rhs) override
    {
      auto& rhs_vec = rhs->getVector();
      auto& u_aux_vec = u_aux->getVector(1);

      for (int i=0; i < rhs_vec.shape()[0]; ++i)
        if (m_block_dependent == BlockDependent::FirstBlock)
        {
          rhs_vec[i] = 2*t + 0.2*std::sqrt(u_aux_vec[i]);
        } else
        {
          rhs_vec[i] = 2*t;
        }

      rhs->markVectorModified();
    }

    void computeJacobian(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, const Real t, linear_system::AssemblerPtr assembler) override
    {
      // nothing to do: d rhs/du = 0
    }

    void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out) override
    {
      auto& v_in  = vec_in->getVector();
      auto& v_out = vec_out->getVector();
      for (int i=0; i < v_in.shape()[0]; ++i)
        v_out[i] = v_in[i];
    }

    void computeMassMatrix(linear_system::AssemblerPtr assembler) override
    {
      auto vol_disc            = getDiscretization()->getVolDisc(0);
      auto dof_numbering       = getDiscretization()->getDofNumbering();
      auto& dof_nums           = dof_numbering->getDofs(0);
      int num_dofs_per_element = vol_disc->getNumSolPtsPerElement();
      int num_el               = vol_disc->getNumElems();

      std::vector<DofInt> dofs(num_dofs_per_element);
      ArrayType<Real, 2> jac(boost::extents[num_dofs_per_element][num_dofs_per_element]);
      std::set<DofInt> seen_dofs;
      for (int i=0; i < num_el; ++i)
      {
        for (int j=0; j < num_dofs_per_element; ++j)
        {
          dofs[j] = dof_nums[i][j];
          jac[j][j] = seen_dofs.count(dofs[j]) == 0 ? 1 : 0;
          seen_dofs.insert(dofs[j]);
        }

        assembler->assembleVolume(0, i, jac);
      }
    }

    AuxiliaryEquationsPtr getAuxEquations() override { return m_aux_eqns; }

  private:
    AuxiliaryEquationsPtr m_aux_eqns;
    BlockDependent m_block_dependent;
};

}

TEST_F(CNDependentTester2, SecondBlockDependentExactnessness)
{

  // The exact solution is u(t) = t^2, and u_aux(t) = t^2 + 0.1*t^2 + 1
  // Note that d u_aux/dt  = 2*t + 0.2*t = 2*t + 0.2*sqrt(u(t)), so this
  // test verifies that the correct u vector is being passed throughout the 
  // Newton/Crank-Nicolson functions
  Real delta_t = 0.5;
  timesolvers::TimeStepperOpts opts;
  opts.t_start = 1.0;
  opts.t_end   = opts.t_start + 2.5*delta_t;
  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerConstant>(delta_t);
  opts.mat_type = linear_system::LargeMatrixType::Dense;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOpts>();
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 5;  //TODO: test 1

  auto model = std::make_shared<PhysicsModelTest>(disc, BlockDependent::SecondBlock);
  auto u = makeDiscVector(disc);
  auto u_aux = makeAuxiliaryEquationsStorage(model->getAuxEquations());

  u->set(1);
  auto& u_aux_vec = u_aux->getVector(1);
  for (int i=0; i < u_aux_vec.shape()[0]; ++i)
    u_aux_vec[i] = 2.1;

  timesolvers::CrankNicolson crank(model, u, u_aux, opts);
  crank.solve();

  auto& u_vec = u->getVector();
  for (int i=0; i < u_vec.shape()[0]; ++i)
    EXPECT_NEAR(u_vec[i], opts.t_end * opts.t_end, 1e-11);

  for (int i=0; i < u_aux_vec.shape()[0]; ++i)
    u_aux_vec[i] = 1.1*opts.t_end*opts.t_end + 1;
}


TEST_F(CNDependentTester2, FirstBlockDependentExactnessness)
{

  // Similiar to the SecondBlockDependentExactness, except that
  // u is a function of t and u_aux and
  // u_aux is a function of t only (the equations are flipped),
  Real delta_t = 0.5;
  timesolvers::TimeStepperOpts opts;
  opts.t_start = 1.0;
  opts.t_end   = opts.t_start + 2.5*delta_t;
  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerConstant>(delta_t);
  opts.mat_type = linear_system::LargeMatrixType::Dense;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOpts>();
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 5;  //TODO: test 1

  auto model = std::make_shared<PhysicsModelTest>(disc, BlockDependent::FirstBlock);
  auto u = makeDiscVector(disc);
  auto u_aux = makeAuxiliaryEquationsStorage(model->getAuxEquations());

  u->set(2.1);
  auto& u_aux_vec = u_aux->getVector(1);
  for (int i=0; i < u_aux_vec.shape()[0]; ++i)
    u_aux_vec[i] = 1;

  timesolvers::CrankNicolson crank(model, u, u_aux, opts);
  crank.solve();

  auto& u_vec = u->getVector();
  for (int i=0; i < u_vec.shape()[0]; ++i)
    EXPECT_NEAR(u_vec[i], 1.1*opts.t_end*opts.t_end + 1, 1e-11);

  for (int i=0; i < u_aux_vec.shape()[0]; ++i)
    u_aux_vec[i] = opts.t_end * opts.t_end;
}