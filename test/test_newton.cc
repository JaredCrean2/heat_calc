#include "gtest/gtest.h"
#include "discretization/disc_vector.h"
#include "discretization/discretization.h"
#include "linear_system/assembler.h"
#include "physics/AuxiliaryEquations.h"
#include "time_solver/newton.h"
#include "mesh_helper.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern_mesh.h"

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

    int getNumBlocks() const override { return m_aux_eqns->getNumBlocks(); }

    // returns the number of variables in the given block
    int getBlockSize(int block) const override { return m_aux_eqns->getBlockSize(block); }

    void computeRhs(int block, DiscVectorPtr u_vec, ArrayType<Real, 1>& rhs) override { return m_aux_eqns->computeRhs(block, u_vec, 0.0, rhs); }

    void computeJacobian(int block, DiscVectorPtr u_vec, linear_system::LargeMatrixPtr mat) override
    {
      auto assembler = std::make_shared<linear_system::SimpleAssembler>(mat);
      return m_aux_eqns->computeJacobian(block, u_vec, 0.0, assembler);
    }

    void multiplyOffDiagonal(int iblock, int jblock, DiscVectorPtr u_vec, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      m_aux_eqns->multiplyOffDiagonal(iblock, jblock, u_vec, 0.0, x, b);
    }

    void setBlockSolution(int block, const ArrayType<Real, 1>& vals) override
    {
      m_aux_eqns->setBlockSolution(block, vals);
    }

    ArrayType<Real, 1>& getBlockSolution(int block) override
    {
      return m_aux_eqns->getBlockSolution(block);
    }

    AuxiliaryEquationsJacobiansPtr getJacobians() override
    {
      return m_aux_eqns->getJacobians();
    }

    AuxiliaryEquationsStoragePtr createStorage() override
    {
      return std::make_shared<AuxiliaryEquationStorage>(*m_aux_eqns);
    }

  private:
    AuxiliaryEquationsPtr m_aux_eqns;
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

    void resetForNewSolve() override {}
  
    Real computeFunc(const DiscVectorPtr u, bool compute_norm, DiscVectorPtr f) override
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

    void computeJacobian(const DiscVectorPtr u, linear_system::LargeMatrixPtr jac) override
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
  auto result = newton.solve(u, abstol, -1, itermax);

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