#include "gtest/gtest.h"
#include "mesh_helper.h"
#include "physics/PhysicsModel.h"
#include "time_solver/crank_nicolson.h"
#include "physics/heat/HeatEquation.h"
#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC_defs.h"

namespace {


  Real ex_sol(Real x, Real y, Real z, Real t)
  { 
    return 1;
  }

  std::array<Real, 3> deriv(Real x, Real y, Real z, Real t)
  { 
    return std::array<Real, 3>{0, 0, 0}; 
  };

  Real src(Real x, Real y, Real z, Real t)
  { 
    return 0; 
  };


class CNPhysicsModel : public PhysicsModel
{
  public:
    explicit CNPhysicsModel(std::shared_ptr<Heat::HeatEquation> heat) :
      PhysicsModel(heat->getDiscretization()),
      m_heat(heat)
    {}

    virtual ~CNPhysicsModel() {}

    // overwrites rhs with the right hand side
    // on entry, u has the solution in vector form
    // on exit, rhs has the residual in array form
    virtual void computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs)
    {
      //m_heat->computeRhs(u, t, rhs);
      rhs->set(0);
      //rhs->syncArrayToVector();
      auto& vec = rhs->getVector();
      for (int i=0; i < vec.shape()[0]; ++i)
        vec[i] += t;
    }

    virtual void computeJacobian(DiscVectorPtr u, const Real t, linear_system::AssemblerPtr assembler)
    {
      //m_heat->computeJacobian(u, t, assembler);
    }

    virtual void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out)
    {
      //m_heat->applyMassMatrix(vec_in, vec_out);
      auto& v_vec_in = vec_in->getVector();
      auto& v_vec_out = vec_out->getVector();
      for (int i=0; i < vec_in->getNumDofs(); ++i)
        v_vec_out[i] = v_vec_in[i];
      vec_out->markVectorModified();
    }

    virtual void computeMassMatrix(linear_system::AssemblerPtr assembler)
    {
      //m_heat->computeMassMatrix(assembler);
      std::vector<DofInt> dofs(1);
      ArrayType<Real, 2> vals(boost::extents[1][1]);
      auto mat = assembler->getMatrix();

      for (int i=0; i < mat->getMLocal(); ++i)
      {
        dofs[0]    = i;
        vals[0][0] = 1 * assembler->getAlpha();
        mat->assembleValues(dofs, vals);
      }
    }

    void addDirichletBC(DirichletBCPtr bc)
    {
       PhysicsModel::addDirichletBC(bc);
       m_heat->addDirichletBC(bc); 
    }

    void addNeumannBC(NeumannBCPtr bc) 
    { 
      PhysicsModel::addNeumannBC(bc);
      m_heat->addNeumannBC(bc);
    }

    const std::vector<DirichletBCPtr>& getDirichletBCs() const {return m_heat->getDirichletBCs(); }

    const std::vector<NeumannBCPtr>& getNeumannBCs() const { return m_heat->getNeumannBCs(); }

    void addSourceTerm(SourceTermPtr src) 
    { 
      PhysicsModel::addSourceTerm(src);
      m_heat->addSourceTerm(src);
    }

    SourceTermPtr getSourceTerm(int idx) const 
    { 
      return m_heat->getSourceTerm(idx);
    }

  private:
    std::shared_ptr<Heat::HeatEquation> m_heat;
};




class CNTester : public StandardDiscSetup,
                 public testing::Test
{
  protected:
    using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

    CNTester()
    {
      Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 2, 2, 2);
      setup(3, 1, spec, {false, false, false, false, false, false});
      setSolution();
    }

    void setSolution(const Heat::VolumeGroupParams& params = Heat::VolumeGroupParams{1, 1, 1})
    {
      auto heat    = std::make_shared<Heat::HeatEquation>(disc);
      cn_model = std::make_shared<CNPhysicsModel>(heat);
      u_vec   = makeDiscVector(disc);
      //res_vec = makeDiscVector(disc);

      for (int i=0; i < disc->getNumVolDiscs(); ++i)
      {
        cn_model->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
        heat->addVolumeGroupParams(params);
      }

      for (int i=0; i < disc->getNumBCSurfDiscs(); ++i)
      {
        auto surf = disc->getBCSurfDisc(i);
        if (surf->getIsDirichlet())
          cn_model->addDirichletBC(makeDirichletBCMMS(surf, ex_sol));
        else
          cn_model->addNeumannBC(makeNeumannBCMMS(surf, deriv));
      }


      //res_vec->set(0);
      auto f = [&](Real x, Real y, Real z)
                  { return ex_sol(x, y, z, 0); };
      u_vec->setFunc(f);
    }

    std::shared_ptr<CNPhysicsModel> cn_model;
    DiscVectorPtr u_vec;
    //DiscVectorPtr res_vec;
};
}


TEST_F(CNTester, Linear)
{
  timesolvers::TimeStepperOpts opts;
  opts.t_start = 0.0;
  opts.t_end   = 0.55;
  opts.delta_t = 0.1;
  opts.mat_type = linear_system::LargeMatrixType::Dense;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOpts>();
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 5;  //TODO: test 1

  timesolvers::CrankNicolson crank(cn_model, u_vec, opts);
  crank.solve();

  if (!u_vec->isVectorCurrent())
    u_vec->syncArrayToVector();

  auto& vec = u_vec->getVector();
  for (int i=0; i < vec.shape()[0]; ++i)
  {
    std::cout << "dof " << i << std::endl;
    EXPECT_NEAR(vec[i], 0.5*opts.t_end*opts.t_end + 1, 1e-12);
  }

}