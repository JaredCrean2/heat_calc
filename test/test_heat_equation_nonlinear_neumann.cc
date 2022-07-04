#include "gtest/gtest.h"
#include "discretization/NeumannBC.h"
#include "discretization/surface_discretization.h"
#include "mesh_helper.h"
#include "test_helper.h"
#include "physics/heat/bc_defs.h"

#include <random>

#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC_defs.h"
#include "discretization/disc_vector.h"
#include "discretization/source_term.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/basis_vals.h"
#include "linear_system/large_matrix_dense.h"

namespace {

class NeumannBCTester : public StandardDiscSetup,
                        public ::testing::Test
{
  public:

    NeumannBCTester()
    {
      setup(3, 1, getStandardMeshSpec(), {false, false, false, false, false, false});
    }

    void testDerivative(NeumannBCPtr bc)
    {
      auto surf = bc->getSurfDisc();
      std::vector<Real> u_vals(surf->getNumQuadPtsPerFace()), flux_vals_deriv(3*surf->getNumQuadPtsPerFace());
      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        u_vals[i] = 1 + 0.1 * i;

      bc->getValueDeriv(0, 0.0, u_vals.data(), flux_vals_deriv.data());
      auto flux_vals_deriv_fd = computeDerivFiniteDifference(bc, u_vals);

      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
        EXPECT_NEAR(flux_vals_deriv[i], flux_vals_deriv_fd[i], 1e-6);
    }

  private:

    std::vector<Real> computeDerivFiniteDifference(NeumannBCPtr bc, std::vector<Real> u_vals)
    {
      auto surf = bc->getSurfDisc();
      std::vector<Real> flux_vals(3*surf->getNumQuadPtsPerFace()), 
                        flux_vals2(3*surf->getNumQuadPtsPerFace()),
                        flux_vals_deriv_fd(3*surf->getNumQuadPtsPerFace());

      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
      {
        flux_vals[i]  = 666 + i;
        flux_vals2[i] = 666 + i;
      }

      bc->getValue(0, 0.0, u_vals.data(), flux_vals.data());
      Real eps = 1e-7;
      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        u_vals[i] += eps;

      bc->getValue(0, 0.0, u_vals.data(), flux_vals2.data());
      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
          flux_vals_deriv_fd[i] = (flux_vals2[i] - flux_vals[i])/eps;

      return flux_vals_deriv_fd;
    }
};


class NeumannBCQuadratic : public NeumannBC
{
  public:
    NeumannBCQuadratic(SurfDiscPtr surf) :
      NeumannBC(surf, true)
    {}

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
        for (int d=0; d < 3; ++d)
          flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = (sol_vals[i] + 0.2 * d) * (sol_vals[i] + 0.2 * d);
    }

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
        for (int d=0; d < 3; ++d)
          flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = 2 * (sol_vals[i] + 0.2 * d);
    }
};
}

//-----------------------------------------------------------------------------
// Test individual Jacobians


TEST_F(NeumannBCTester, QuadraticBC)
{
  auto bc = std::make_shared<NeumannBCQuadratic>(disc->getSurfDisc(0));
  testDerivative(bc);
}

TEST_F(NeumannBCTester, NewtonCoolingBC)
{
  auto bc = std::make_shared<Heat::NewtonCooling>(disc->getSurfDisc(0), 2);
  bc->setExternalTemperature(5);
  testDerivative(bc);
}


//-----------------------------------------------------------------------------
// Test overall Jacobian


namespace {

  Real zeroFunc(Real x, Real y, Real z, Real t)
  {
    return 0;
  }

  class HeatNeumannTester : public StandardDiscSetup,
                            public testing::Test
  {
    protected:
      using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

      HeatNeumannTester()
      {
        setup();
      }

      template <typename Tex>
      void setSolution(Tex ex_sol, const Heat::VolumeGroupParams& params = Heat::VolumeGroupParams{1, 1, 1})
      {
        heat = std::make_shared<Heat::HeatEquation>(disc);
        u_vec = makeDiscVector(disc);
        res_vec = makeDiscVector(disc);


        for (int i=0; i < disc->getNumVolDiscs(); ++i)
        {
          heat->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), zeroFunc));
          heat->addVolumeGroupParams(params);
        }

        for (int i=0; i < disc->getNumBCSurfDiscs(); ++i)
        {
          auto surf = disc->getBCSurfDisc(i);
          heat->addNeumannBC(std::make_shared<NeumannBCQuadratic>(surf));
        }


        res_vec->set(0);
        auto f = [&](Real x, Real y, Real z)
                    { return ex_sol(x, y, z, 0); };
        u_vec->setFunc(f);
      }

      HeatPtr heat;
      DiscVectorPtr u_vec;
      DiscVectorPtr res_vec;
  };
}



TEST_F(HeatNeumannTester, JacobianFiniteDifferenceDirichlet)
{
  SERIAL_ONLY();

  // Note: this tests Dirichlet surfaces and Neumann surfaces where the flux
  //       does not depend on the solution
  using Rng = std::mt19937;
  Rng rng;
  const int seed = 42;
  rng.seed(seed);
  const Real eps = 1e-6;
  std::uniform_real_distribution<Real> uniform_rng(-1, 1);
  const int nvectors = 10;

  linear_system::LargeMatrixOpts opts;
  opts.factor_in_place = false;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric = false;
  std::vector<bool> dirichlet_surfs = {false, false, false, false, false, false};
  int sol_degree = 1;

  auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                      { return x + 2*y * 3*z; };

  setup(2*sol_degree, sol_degree, dirichlet_surfs);

  auto num_dofs  = disc->getDofNumbering()->getNumOwnedDofs();
  auto mat       = std::make_shared<linear_system::LargeMatrixDense>(num_dofs, num_dofs, opts);
  auto assembler = std::make_shared<linear_system::Assembler>(disc, mat);
  auto res_vec2  = makeDiscVector(disc);
  ArrayType<Real, 1> pert_vec(boost::extents[num_dofs]);
  std::vector<Real> product_fd(num_dofs);
  ArrayType<Real, 1>  product_jac(boost::extents[num_dofs]);
  
  setSolution(ex_sol_l);

  heat->computeJacobian(u_vec, 0.0, assembler);


  for (int i=0; i < nvectors; ++i)
  {
    setSolution(ex_sol_l);

    for (int j=0; j < pert_vec.shape()[0]; ++j)
    {
      res_vec->getVector()[j] = 0;
      res_vec2->getVector()[j] = 0;
    }

    // compute at original state
    heat->computeRhs(u_vec, 0.0, res_vec);
    res_vec->syncArrayToVector();

    // apply perturbation
    for (unsigned int j=0; j < pert_vec.shape()[0]; ++j)
      pert_vec[j] = uniform_rng(rng);


    for (int j=0; j < num_dofs; ++j)
      u_vec->getVector()[j] += eps * pert_vec[j];
    u_vec->markVectorModified();

    heat->computeRhs(u_vec, 0.0, res_vec2);
    res_vec2->syncArrayToVector();

    for (int j=0; j < num_dofs; ++j)
      product_fd[j] = (res_vec2->getVector()[j] - res_vec->getVector()[j])/eps;

    // compute explicit jacobian vector product
    mat->matVec(pert_vec, product_jac);

    for (int j=0; j < num_dofs; ++j)
      EXPECT_NEAR(product_fd[j], product_jac[j], 1e-5);
  }
}