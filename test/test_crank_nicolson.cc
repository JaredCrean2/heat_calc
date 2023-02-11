#include "gtest/gtest.h"
#include "discretization/DirichletBC.h"
#include "linear_system/large_matrix_factory.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern_dense.h"
#include "mesh_helper.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/PhysicsModel.h"
#include "test_helper.h"
#include "time_solver/crank_nicolson.h"
#include "physics/heat/HeatEquation.h"
#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC_defs.h"
#include <random>

namespace {


Real ex_sol(Real x, Real y, Real z, Real t, int degree_space, int degree_time)
{
  return std::pow(t, degree_time) + std::pow(x, degree_space) + std::pow(y, degree_space) + std::pow(z, degree_space);
  //return std::pow(y, degree);
}

Real ex_sol_dt(Real x, Real y, Real z, Real t, int degree_space, int degree_time)
{
  Real sol_dt = 0;
  if (degree_time > 0)
    sol_dt = degree_time * std::pow(t, degree_time - 1);

  return sol_dt;
}


std::array<Real, 3> ex_sol_deriv(Real x, Real y, Real z, Real t, int degree_space)
{
  std::array<Real, 3> derivs{0, 0, 0};
  if (degree_space > 0)
  {
    derivs[0] = degree_space * std::pow(x, degree_space - 1);
    derivs[1] = degree_space * std::pow(y, degree_space - 1);
    derivs[2] = degree_space * std::pow(z, degree_space - 1);
  }

  return derivs;
}


Real src_func_dir(Real x, int degree_space)
{
  Real spacial_term = 0;
  if (degree_space >= 2)
    spacial_term = -(degree_space * (degree_space - 1) * std::pow(x, degree_space-2));

  return spacial_term;
}

Real src_func(Real x, Real y, Real z, Real t, int degree_space, int degree_time)
{
  Real src = src_func_dir(x, degree_space) + src_func_dir(y, degree_space) + src_func_dir(z, degree_space);

  if (degree_time >= 1)
    src += degree_time * std::pow(t, degree_time - 1);

  return src;
}

enum SolveType
{
  TIMEONLY,
  SPACETIME
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
    virtual void computeRhs(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, const Real t, DiscVectorPtr rhs)
    {
      if (m_solve_type == SPACETIME)
      {
        m_heat->computeRhs(u, u_aux, t, rhs);
      } else
      {
        rhs->set(0);
        auto& vec = rhs->getVector();
        for (int i=0; i < vec.shape()[0]; ++i)
          vec[i] += t;
      }
    }

    virtual void computeJacobian(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, const Real t, linear_system::AssemblerPtr assembler)
    {
      if (m_solve_type == SPACETIME)
        m_heat->computeJacobian(u, u_aux, t, assembler);
    }

    virtual void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out)
    {
      if (m_solve_type == SPACETIME)
      {
        m_heat->applyMassMatrix(vec_in, vec_out);
      } else
      {
        auto& v_vec_in = vec_in->getVector();
        auto& v_vec_out = vec_out->getVector();
        for (int i=0; i < vec_in->getNumDofs(); ++i)
          v_vec_out[i] = v_vec_in[i];
        vec_out->markVectorModified();
      }
    }

    virtual void computeMassMatrix(linear_system::AssemblerPtr assembler)
    {
      if (m_solve_type == SPACETIME)
      {
        m_heat->computeMassMatrix(assembler);
      } else
      {
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

    void setSolveType(SolveType solve_type) { m_solve_type = solve_type;}

  private:
    std::shared_ptr<Heat::HeatEquation> m_heat;
    SolveType m_solve_type = TIMEONLY;
};




class CNTester : public StandardDiscSetup,
                 public testing::Test
{
  protected:
    using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

    CNTester()
    {
      Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 3, 3, 3);
      setup(3, 1, spec, {false, false, false, false, false, false});
      //setup(3, 1, spec, {true, true, true, true, true, true});

    }

    template <typename Tex, typename Tderiv, typename Tsrc, typename Tsrc2 = impl::ZeroFuncType>
    void setSolution(Tex ex_sol, Tderiv deriv, Tsrc src, Tsrc2 ex_sol_dt=&(impl::zeroFunc), const Heat::VolumeGroupParams& params = Heat::VolumeGroupParams{1, 1, 1})
    {
      auto heat    = std::make_shared<Heat::HeatEquation>(disc);
      cn_model = std::make_shared<CNPhysicsModel>(heat);
      u_vec   = makeDiscVector(disc);
      u_aux_vec = makeAuxiliaryEquationsStorage(heat->getAuxEquations());

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
          cn_model->addDirichletBC(makeDirichletBCMMS(surf, ex_sol, ex_sol_dt));
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
    AuxiliaryEquationsStoragePtr u_aux_vec;
    //DiscVectorPtr res_vec;
};

linear_system::LargeMatrixOptsPetsc get_options()
{
  linear_system::LargeMatrixOptsPetsc opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;
  opts.petsc_opts["ksp_atol"] = "1e-13";
  opts.petsc_opts["ksp_rtol"] = "1e-50";
  opts.petsc_opts["ksp_monitor"] = "";

  return opts;
}

}


TEST_F(CNTester, Linear)
{
  SERIAL_ONLY();
  timesolvers::TimeStepperOpts opts;
  opts.t_start = 0.0;
  opts.t_end   = 0.55;
  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerConstant>(0.1);
  opts.mat_type = linear_system::LargeMatrixType::Dense;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOpts>();
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 5;  //TODO: test 1

  Real kappa = 1;
  int degree = 0;
  auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                      { return ex_sol(x, y, z, t, degree, degree); };

  auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> std::array<Real, 3>
                      { 
                        auto vals = ex_sol_deriv(x, y, z, t, degree);
                        for (auto& v : vals)
                          v *= kappa;
                        return vals;
                      };
  auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { return kappa * src_func(x, y, z, t, degree, degree); };

  setSolution(ex_sol_l, deriv_l, src_func_l);

  std::cout << "u_aux_vec = " << u_aux_vec << std::endl;
  timesolvers::CrankNicolson crank(cn_model, u_vec, u_aux_vec, opts);
  crank.solve();

  if (!u_vec->isVectorCurrent())
    u_vec->syncArrayToVector();

  auto& vec = u_vec->getVector();
  for (int i=0; i < vec.shape()[0]; ++i)
  {
    EXPECT_NEAR(vec[i], 0.5*opts.t_end*opts.t_end + 4, 1e-12);
  }
}

// This test has to use Neumann BCs because we don't compute
// the initial condition the way the theory says we should
// for time-varying Dirichlet BCs

TEST_F(CNTester, PolynomialExactness)
{
  timesolvers::TimeStepperOpts opts;
  opts.t_start = 0.0;
  opts.t_end   = 0.55;
  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerConstant>(0.025);
  opts.mat_type = linear_system::LargeMatrixType::Petsc;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>(get_options());
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 1;

  Real kappa = 1;
  int degree_time = 2;

  //Not sure why this is exact up to degree 4
  for (int degree_space=0; degree_space < 5; ++degree_space)
  {
    std::cout << "testing degree_space = " << degree_space << std::endl;
    auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { return ex_sol(x, y, z, t, degree_space, degree_time); };

    auto ex_sol_dt_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { return ex_sol_dt(x, y, z, t, degree_space, degree_time); };

    auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> std::array<Real, 3>
                        { 
                          auto vals = ex_sol_deriv(x, y, z, t, degree_space);
                          for (auto& v : vals)
                            v *= kappa;
                          return vals;
                        };
    auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return kappa * src_func(x, y, z, t, degree_space, degree_time); };

    setSolution(ex_sol_l, deriv_l, src_func_l, ex_sol_dt_l);
    cn_model->setSolveType(SPACETIME);

    timesolvers::CrankNicolson crank(cn_model, u_vec, u_aux_vec, opts);
    crank.solve();

    if (!u_vec->isArrayCurrent())
      u_vec->syncVectorToArray();

    auto disc = cn_model->getDiscretization();
    for (int i=0; i < disc->getNumVolDiscs(); ++i)
    {
      auto vol_disc = disc->getVolDisc(i);
      auto& coords  = vol_disc->vol_group.coords;
      auto& u_arr   = u_vec->getArray(i);
      for (int el=0; el < vol_disc->getNumElems(); ++el)
        for (int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
        {
          Real ex_val = ex_sol_l(coords[el][j][0], coords[el][j][1], coords[el][j][2], opts.t_end);
          //std::cout << "ex_val = " << ex_val << ", u_val = " << u_arr[el][j] << std::endl;
          EXPECT_NEAR(u_arr[el][j], ex_val, 1e-12);
        }
    }
  }
}

//TODO: write test verifying consistency of Jacobian and residual
TEST_F(CNTester, JacobianFD)
{
  SERIAL_ONLY();

  Real kappa = 1;
  int degree_space = 2, degree_time = 2;

  auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                      { return ex_sol(x, y, z, t, degree_space, degree_time); };

  auto ex_sol_dt_l = [&] (Real x, Real y, Real z, Real t) -> Real
                      { return ex_sol_dt(x, y, z, t, degree_space, degree_time); };

  auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> std::array<Real, 3>
                      { 
                        auto vals = ex_sol_deriv(x, y, z, t, degree_space);
                        for (auto& v : vals)
                          v *= kappa;
                        return vals;
                      };
  auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { return kappa * src_func(x, y, z, t, degree_space, degree_time); };
  setSolution(ex_sol_l, deriv_l, src_func_l, ex_sol_dt_l);

  cn_model->setSolveType(SPACETIME);
  const auto num_dofs  = disc->getDofNumbering()->getNumOwnedDofs();
  auto mat_opts = std::make_shared<linear_system::LargeMatrixOpts>();
  auto mat = linear_system::largeMatrixFactory(linear_system::LargeMatrixType::Dense,
                                               mat_opts, std::make_shared<linear_system::SparsityPatternDense>(num_dofs));
  const Real tn = 1.0;
  const Real delta_t = 0.5;
  timesolvers::CrankNicolsonFunction cn_func(cn_model, mat, tn - delta_t);

  using Rng = std::mt19937;
  Rng rng;
  const int seed = 42;
  const Real eps = 1e-6;
  std::uniform_real_distribution<Real> uniform_rng(-1, 1);
  const int nvectors = 10;
  ArrayType<Real, 1> pert_vec(boost::extents[num_dofs]);
  std::vector<Real> product_fd(num_dofs);
  ArrayType<Real, 1>  product_jac(boost::extents[num_dofs]);
  rng.seed(seed);
  auto res_vec1 = makeDiscVector(disc);
  auto res_vec2 = makeDiscVector(disc);


  // set previous timestep solution
  Real t = tn - delta_t;
  auto ex_sol_t = [&] (Real x, Real y, Real z) -> Real
                      { return ex_sol_l(x, y, z, t); };

  u_vec->setFunc(ex_sol_t);
  cn_func.setTnp1(u_vec, u_aux_vec, tn);

  // get current timestep residual and jacobian
  t = tn;
  u_vec->setFunc(ex_sol_t);
  cn_func.computeFunc(u_vec, u_aux_vec, t,res_vec1);
  cn_func.computeJacobian(u_vec, u_aux_vec, mat);

  for (int i=0; i < nvectors; ++i)
  {
    u_vec->setFunc(ex_sol_t);

    // apply perturbation
    for (unsigned int j=0; j < pert_vec.shape()[0]; ++j)
    {
      pert_vec[j] = uniform_rng(rng);
    }

    for (int j=0; j < num_dofs; ++j)
      u_vec->getVector()[j] += eps * pert_vec[j];
    u_vec->markVectorModified();

    cn_func.computeFunc(u_vec, u_aux_vec, false, res_vec2);

    for (int j=0; j < num_dofs; ++j)
      product_fd[j] = (res_vec2->getVector()[j] - res_vec1->getVector()[j])/eps;

    mat->matVec(pert_vec, product_jac);

    for (int j=0; j < num_dofs; ++j)
      EXPECT_NEAR(product_fd[j], product_jac[j], 1e-5);
  }
}