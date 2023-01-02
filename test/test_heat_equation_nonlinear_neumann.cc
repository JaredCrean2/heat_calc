#include "gtest/gtest.h"
#include "discretization/NeumannBC.h"
#include "discretization/surface_discretization.h"
#include "mesh_helper.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/heat/solar_position.h"
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

    void testDerivativeWrtTair(std::shared_ptr<Heat::AirWindSkyNeumannBC> bc)
    {
      Real air_temp = 300;
      Real air_speed = 5;
      std::array<Real, 3> air_direction{1.0/std::sqrt(14.0), 2.0/std::sqrt(14.0), 3.0/std::sqrt(14.0)};
      Real ir_horizontal_radiation = 10;
      Real direct_normal_radiation = 20;
      Real diffuse_radiation = 30;
      Heat::DirectionCosines solar_direction{std::cos(0.25), std::cos(0.5), std::cos(0.9)};

      bc->setAirTemperature(air_temp);
      bc->setAirSpeed(air_speed);
      bc->setAirDirection(air_direction);
      bc->setIRHorizontalRadiation(ir_horizontal_radiation);
      bc->setDirectNormalRadiation(direct_normal_radiation);
      bc->setDiffuseRadiation(diffuse_radiation);
      bc->setSolarDirection(solar_direction);

      auto surf = bc->getSurfDisc();
      std::vector<Real> u_vals(surf->getNumQuadPtsPerFace()), flux_vals_deriv(3*surf->getNumQuadPtsPerFace());
      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        u_vals[i] = 1 + 0.1 * i;

      auto dflux_dtair_fd = computeTairFiniteDifference(bc, u_vals, air_temp);

      std::vector<Real> flux(3*surf->getNumQuadPtsPerFace()), dflux_dtair(3*surf->getNumQuadPtsPerFace());
      bc->setAirTemperature(air_temp);
      bc->getValuedTair(0, 0.0, u_vals.data(), flux.data(), dflux_dtair.data());
      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
        EXPECT_NEAR(dflux_dtair[i], dflux_dtair_fd[i], 5e-5);
    }

    void testDerivativeWrtTwall_rev(std::shared_ptr<Heat::AirWindSkyNeumannBC> bc)
    {
      //TODO: move setup to own function
      Real air_temp = 300;
      Real air_speed = 5;
      std::array<Real, 3> air_direction{1.0/std::sqrt(14.0), 2.0/std::sqrt(14.0), 3.0/std::sqrt(14.0)};
      Real ir_horizontal_radiation = 10;
      Real direct_normal_radiation = 20;
      Real diffuse_radiation = 30;
      Heat::DirectionCosines solar_direction{std::cos(0.25), std::cos(0.5), std::cos(0.9)};

      bc->setAirTemperature(air_temp);
      bc->setAirSpeed(air_speed);
      bc->setAirDirection(air_direction);
      bc->setIRHorizontalRadiation(ir_horizontal_radiation);
      bc->setDirectNormalRadiation(direct_normal_radiation);
      bc->setDiffuseRadiation(diffuse_radiation);
      bc->setSolarDirection(solar_direction);

      // test against forward mode code
      auto surf = bc->getSurfDisc();
      std::vector<Real> u_vals(surf->getNumQuadPtsPerFace()), u_vals_bar(surf->getNumQuadPtsPerFace()),
                        flux_vals_deriv(3*surf->getNumQuadPtsPerFace()),
                        flux_vals_deriv_rev(3*surf->getNumQuadPtsPerFace()),
                        flux_vals_bar(3*surf->getNumQuadPtsPerFace());

      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        u_vals[i] = 1 + 0.1 * i;

      bc->getValueDeriv(0, 0.0, u_vals.data(), flux_vals_deriv.data());


      for (int d=0; d < 3; ++d)
      {
        std::fill(u_vals_bar.begin(), u_vals_bar.end(), 0);
        std::fill(flux_vals_bar.begin(), flux_vals_bar.end(), 0);
        for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
          flux_vals_bar[d * surf->getNumQuadPtsPerFace() + i] = 1;

        bc->getValue_rev(0, 0.0, u_vals.data(), u_vals_bar.data(), flux_vals_bar.data());

        for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
          flux_vals_deriv_rev[d * surf->getNumQuadPtsPerFace() + i] = u_vals_bar[i];
      }

      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
        EXPECT_NEAR(flux_vals_deriv[i], flux_vals_deriv_rev[i], 1e-13);
      
    }

  private:

    void setParameters(std::shared_ptr<Heat::AirWindSkyNeumannBC> bc)
    {
      Real air_temp = 300;
      Real air_speed = 5;
      std::array<Real, 3> air_direction{1.0/std::sqrt(14.0), 2.0/std::sqrt(14.0), 3.0/std::sqrt(14.0)};
      Real ir_horizontal_radiation = 10;
      Real direct_normal_radiation = 20;
      Real diffuse_radiation = 30;
      Heat::DirectionCosines solar_direction{std::cos(0.25), std::cos(0.5), std::cos(0.9)};

      bc->setAirTemperature(air_temp);
      bc->setAirSpeed(air_speed);
      bc->setAirDirection(air_direction);
      bc->setIRHorizontalRadiation(ir_horizontal_radiation);
      bc->setDirectNormalRadiation(direct_normal_radiation);
      bc->setDiffuseRadiation(diffuse_radiation);
      bc->setSolarDirection(solar_direction);
    }

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

    std::vector<Real> computeTairFiniteDifference(std::shared_ptr<Heat::AirWindSkyNeumannBC> bc, std::vector<Real> u_vals, Real t_air)
    {
      auto surf = bc->getSurfDisc();
      std::vector<Real> flux_vals(3*surf->getNumQuadPtsPerFace()), 
                        flux_vals2(3*surf->getNumQuadPtsPerFace()),
                        flux_vals_deriv_fd(3*surf->getNumQuadPtsPerFace());

      Real eps = 1e-7;
      bc->setAirTemperature(t_air);
      bc->getValue(0, 0.0, u_vals.data(), flux_vals.data());

      bc->setAirTemperature(t_air + eps);
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
  SERIAL_ONLY()
  auto bc = std::make_shared<NeumannBCQuadratic>(disc->getSurfDisc(0));
  testDerivative(bc);
}

TEST_F(NeumannBCTester, NewtonCoolingBC)
{
  SERIAL_ONLY()
  auto bc = std::make_shared<Heat::NewtonCooling>(disc->getSurfDisc(0), 2);
  bc->setExternalTemperature(5);
  testDerivative(bc);
}

TEST_F(NeumannBCTester, TarpBC)
{
  SERIAL_ONLY()
  Real air_temp       = 20;    
  Real surface_area   = 4;
  Real perimeter      = 8;
  std::array<Real, 3> vertical_vector{0, 1, 0};
  std::array<Real, 3> pt_at_zero_altitude{-1, -1, -1};
  std::array<Real, 3> air_direction{0, 1, 0};
  int met_terrain_index   = 0;
  Real met_altitude       = 10;
  int local_terrain_index = 1;
  Real air_speed          = 0;
  int roughness_type = 0;

  auto bc = std::make_shared<Heat::TarpBC>(disc->getSurfDisc(0), surface_area, perimeter, 
                 roughness_type, vertical_vector, pt_at_zero_altitude, met_terrain_index,
                                           met_altitude, local_terrain_index);

  bc->setAirTemperature(air_temp);
  bc->setAirSpeed(air_speed);
  bc->setAirDirection(air_direction);

  testDerivative(bc);
  testDerivativeWrtTair(bc);
  testDerivativeWrtTwall_rev(bc);
}


TEST_F(NeumannBCTester, SkyRadiationBC)
{
  SERIAL_ONLY()
  Real emissivity = 0.9;
  Real sigma = 5.6697e-8;
  std::array<Real, 3> vertical{0, 0, 1};
  Real t_air = 200;
  Real t_sky = 300;
  Real ir_radiation = sigma * std::pow(t_sky, 4);

  auto bc = std::make_shared<Heat::SkyRadiationBC>(disc->getSurfDisc(0), emissivity, vertical);
  bc->setIRHorizontalRadiation(ir_radiation);
  bc->setAirTemperature(t_air);

  testDerivative(bc);
  testDerivativeWrtTair(bc);
  testDerivativeWrtTwall_rev(bc);
}

TEST_F(NeumannBCTester, SolarRadiationBC)
{
  Real absorbtivity = 0.9;
  auto bc = std::make_shared<Heat::SolarRadiationBC>(disc->getSurfDisc(0), absorbtivity);

  Real direct_normal_radiation = 100;
  Real diffuse_radiation = 200;
  Heat::DirectionCosines solar_direction{std::cos(0.25), std::cos(0.5), std::cos(0.9)};

  bc->setDirectNormalRadiation(direct_normal_radiation);
  bc->setDiffuseRadiation(diffuse_radiation);
  bc->setSolarDirection(solar_direction);

  testDerivative(bc);
  testDerivativeWrtTair(bc);
  testDerivativeWrtTwall_rev(bc);
}

TEST_F(NeumannBCTester, SimpleConvectionBC)
{
  Real heat_transfer_coeff = 2;
  auto bc = std::make_shared<Heat::SimpleConvectionBC>(disc->getSurfDisc(0), heat_transfer_coeff);

  testDerivative(bc);
  testDerivativeWrtTair(bc);
  testDerivativeWrtTwall_rev(bc);
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
        u_aux_vec = makeAuxiliaryEquationsStorage(heat->getAuxEquations());
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
      AuxiliaryEquationsStoragePtr u_aux_vec;
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

  heat->computeJacobian(u_vec, u_aux_vec, 0.0, assembler);


  for (int i=0; i < nvectors; ++i)
  {
    setSolution(ex_sol_l);

    for (int j=0; j < pert_vec.shape()[0]; ++j)
    {
      res_vec->getVector()[j] = 0;
      res_vec2->getVector()[j] = 0;
    }

    // compute at original state
    heat->computeRhs(u_vec, u_aux_vec, 0.0, res_vec);
    res_vec->syncArrayToVector();

    // apply perturbation
    for (unsigned int j=0; j < pert_vec.shape()[0]; ++j)
      pert_vec[j] = uniform_rng(rng);


    for (int j=0; j < num_dofs; ++j)
      u_vec->getVector()[j] += eps * pert_vec[j];
    u_vec->markVectorModified();

    heat->computeRhs(u_vec, u_aux_vec, 0.0, res_vec2);
    res_vec2->syncArrayToVector();

    for (int j=0; j < num_dofs; ++j)
      product_fd[j] = (res_vec2->getVector()[j] - res_vec->getVector()[j])/eps;

    // compute explicit jacobian vector product
    mat->matVec(pert_vec, product_jac);

    for (int j=0; j < num_dofs; ++j)
      EXPECT_NEAR(product_fd[j], product_jac[j], 1e-5);
  }
}