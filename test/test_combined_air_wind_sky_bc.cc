#include "gtest/gtest.h"
#include "discretization/surface_discretization.h"
#include "mesh_helper.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/solar_position.h"

namespace {

class AirWindSkyTestBC : public Heat::AirWindSkyNeumannBC
{
  public:
    AirWindSkyTestBC(SurfDiscPtr surf, bool is_nonlinear, Real h1, Real h2, Real h3, Real h4, Real h5) :
      AirWindSkyNeumannBC(surf, is_nonlinear),
      m_h1(h1),
      m_h2(h2),
      m_h3(h3),
      m_h4(h4),
      m_h5(h5)
    {}

    void setAirTemperature(Real temp) override { m_air_temp = temp;}

    Real getAirTemperature() const override { return m_air_temp; }
    
    void setAirSpeed(Real velocity) override { m_air_speed = velocity; }
    
    void setAirDirection(std::array<Real, 3> direction) override {m_air_direction = direction; }

    void setIRHorizontalRadiation(Real flux) override { m_ir_horizontal_radiation = flux; }

    void setDirectNormalRadiation(Real flux) override { m_direct_normal_radiation = flux; }

    void setDiffuseRadiation(Real flux) override { m_diffuse_radiation = flux;}

    void setSolarDirection(const Heat::DirectionCosines& cosines) override { m_direction_cosines = cosines; }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        flux_vals[3*i] = m_h1 * (sol_vals[i] - m_air_temp) +
                         m_h2 * (m_air_direction[0] + m_air_direction[1] + m_air_direction[2]) * m_air_speed + 
                         m_h3 * m_ir_horizontal_radiation +
                         m_h4 * m_direct_normal_radiation +
                         m_h5 * m_diffuse_radiation;
        flux_vals[3*i + 1] = flux_vals[3*i];
        flux_vals[3*i + 2] = flux_vals[3*i];
      }
    }

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals_deriv) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        flux_vals_deriv[3*i]     = m_h1;
        flux_vals_deriv[3*i + 1] = m_h1;
        flux_vals_deriv[3*i + 2] = m_h1;


      }      
    }

    // compute derivative of flux_vals wrt air temperature
    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        flux_vals[3*i] = m_h1 * (sol_vals[i] - m_air_temp) +
                         m_h2 * (m_air_direction[0] + m_air_direction[1] + m_air_direction[2]) * m_air_speed + 
                         m_h3 * m_ir_horizontal_radiation +
                         m_h4 * m_direct_normal_radiation +
                         m_h5 * m_diffuse_radiation;
        flux_vals[3*i + 1] = flux_vals[3*i];
        flux_vals[3*i + 2] = flux_vals[3*i];  

        flux_vals_deriv[3*i]     = -m_h1;
        flux_vals_deriv[3*i + 1] = -m_h1;
        flux_vals_deriv[3*i + 2] = -m_h1;
      }      
    }

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        for (int j=0; j < 3; ++j)
          sol_vals_bar[i] = m_h1 * flux_vals_bar[3*i + j];
      }      
    }

  private:
    Real m_h1, m_h2, m_h3, m_h4, m_h5;
    Real m_air_temp;
    Real m_air_speed;
    std::array<Real, 3> m_air_direction;
    Real m_ir_horizontal_radiation;
    Real m_direct_normal_radiation;
    Real m_diffuse_radiation;
    Heat::DirectionCosines m_direction_cosines;
};

class CombinedAirWindSkyBCTester : public StandardDiscSetup,
                           public testing::Test
{
  protected:
    CombinedAirWindSkyBCTester()
    {
      Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 3, 3, 3);
      setup(3, 1, spec, {false, false, false, false, false, false});
      //setup(3, 1, spec, {true, true, true, true, true, true});
    }

    void setupBCs(const std::vector<Real>& h_vals_first, const std::vector<Real>& h_vals_second)
    {
      surf = disc->getSurfDisc(0);
      m_bc1 = std::make_shared<AirWindSkyTestBC>(surf, false, 
                  h_vals_first[0], h_vals_first[1], h_vals_first[2], h_vals_first[3], h_vals_first[4]);
      m_bc2 = std::make_shared<AirWindSkyTestBC>(surf, false, 
                  h_vals_second[0], h_vals_second[1], h_vals_second[2], h_vals_second[3], h_vals_second[4]);
      m_bc3 = std::make_shared<Heat::CombinedAirWindSkyNeumannBC>(std::vector<std::shared_ptr<Heat::AirWindSkyNeumannBC>>{m_bc1, m_bc2});
    }

    void checkFlux()
    {
      std::vector<Real> sol_vals(surf->getNumQuadPtsPerFace()), 
                        flux_vals1(3*surf->getNumQuadPtsPerFace()),
                        flux_vals2(3*surf->getNumQuadPtsPerFace()),
                        flux_vals3(3*surf->getNumQuadPtsPerFace());

      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        sol_vals[i] = i + 1;

      m_bc1->getValue(0, 0.0, sol_vals.data(), flux_vals1.data());
      m_bc2->getValue(0, 0.0, sol_vals.data(), flux_vals2.data());
      m_bc3->getValue(0, 0.0, sol_vals.data(), flux_vals3.data());

      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
        EXPECT_NEAR(flux_vals3[i], flux_vals1[i] + flux_vals2[i], 1e-13);
    }

    void checkFluxDeriv()
    {
      std::vector<Real> sol_vals(surf->getNumQuadPtsPerFace()), 
                        flux_vals1_deriv(3*surf->getNumQuadPtsPerFace()),
                        flux_vals2_deriv(3*surf->getNumQuadPtsPerFace()),
                        flux_vals3_deriv(3*surf->getNumQuadPtsPerFace());

      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        sol_vals[i] = i + 1;

      m_bc1->getValueDeriv(0, 0.0, sol_vals.data(), flux_vals1_deriv.data());
      m_bc2->getValueDeriv(0, 0.0, sol_vals.data(), flux_vals2_deriv.data());
      m_bc3->getValueDeriv(0, 0.0, sol_vals.data(), flux_vals3_deriv.data());

      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
        EXPECT_NEAR(flux_vals3_deriv[i], flux_vals1_deriv[i] + flux_vals2_deriv[i], 1e-13);
    }   

    void checkFluxDerivTair()
    {
      std::vector<Real> sol_vals(surf->getNumQuadPtsPerFace()), 
                        flux_vals1(3*surf->getNumQuadPtsPerFace()), flux_vals1_deriv(3*surf->getNumQuadPtsPerFace()),
                        flux_vals2(3*surf->getNumQuadPtsPerFace()), flux_vals2_deriv(3*surf->getNumQuadPtsPerFace()),
                        flux_vals3(3*surf->getNumQuadPtsPerFace()), flux_vals3_deriv(3*surf->getNumQuadPtsPerFace());

      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        sol_vals[i] = i + 1;

      m_bc1->getValuedTair(0, 0.0, sol_vals.data(), flux_vals1.data(), flux_vals1_deriv.data());
      m_bc2->getValuedTair(0, 0.0, sol_vals.data(), flux_vals2.data(), flux_vals2_deriv.data());
      m_bc3->getValuedTair(0, 0.0, sol_vals.data(), flux_vals3.data(), flux_vals3_deriv.data());

      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
      {
        EXPECT_NEAR(flux_vals3[i], flux_vals1[i] + flux_vals2[i], 1e-13);
        EXPECT_NEAR(flux_vals3_deriv[i], flux_vals1_deriv[i] + flux_vals2_deriv[i], 1e-13);
      }
    }     

    void checkFlux_rev()
    {
      std::vector<Real> sol_vals(surf->getNumQuadPtsPerFace()), 
                        flux_vals_bar(3*surf->getNumQuadPtsPerFace()),
                        sol_vals_bar1(surf->getNumQuadPtsPerFace()),
                        sol_vals_bar2(surf->getNumQuadPtsPerFace()),
                        sol_vals_bar3(surf->getNumQuadPtsPerFace());

      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        sol_vals[i]      = i + 1;

      for (int i=0; i < 3*surf->getNumQuadPtsPerFace(); ++i)
        flux_vals_bar[i] = i + 1;


      m_bc1->getValue_rev(0, 0.0, sol_vals.data(), sol_vals_bar1.data(), flux_vals_bar.data());
      m_bc2->getValue_rev(0, 0.0, sol_vals.data(), sol_vals_bar2.data(), flux_vals_bar.data());
      m_bc3->getValue_rev(0, 0.0, sol_vals.data(), sol_vals_bar3.data(), flux_vals_bar.data());

      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
      {
        EXPECT_NEAR(sol_vals_bar3[i], sol_vals_bar1[i] + sol_vals_bar2[i], 1e-13);
      }
    }

    SurfDiscPtr surf;
    std::shared_ptr<Heat::AirWindSkyNeumannBC> m_bc1;
    std::shared_ptr<Heat::AirWindSkyNeumannBC> m_bc2;
    std::shared_ptr<Heat::CombinedAirWindSkyNeumannBC> m_bc3;
};
}

TEST_F(CombinedAirWindSkyBCTester, Delta_T)
{
  setupBCs({1, 0, 0, 0, 0}, {2, 0, 0, 0, 0});
  checkFlux();
  checkFluxDeriv();
  checkFluxDerivTair();
  checkFlux_rev();
}

TEST_F(CombinedAirWindSkyBCTester, AirSpeed)
{
  setupBCs({0, 1, 0, 0, 0}, {0, 2, 0, 0, 0});
  checkFlux();
  checkFluxDeriv();
  checkFluxDerivTair();
  checkFlux_rev();
}


TEST_F(CombinedAirWindSkyBCTester, IRHorizontalRadiation)
{
  setupBCs({0, 0, 1, 0, 0}, {0, 0, 2, 0, 0});
  checkFlux();
  checkFluxDeriv();
  checkFluxDerivTair();
  checkFlux_rev();
}

TEST_F(CombinedAirWindSkyBCTester, DirectNormalRadiation)
{
  setupBCs({0, 0, 0, 1, 0}, {0, 0, 0, 2, 0});
  checkFlux();
  checkFluxDeriv();
  checkFluxDerivTair();
  checkFlux_rev();
}

TEST_F(CombinedAirWindSkyBCTester, DiffuseRadiation)
{
  setupBCs({0, 0, 0, 0, 1}, {0, 0, 0, 0, 2});
  checkFlux();
  checkFluxDeriv();
  checkFluxDerivTair();
  checkFlux_rev();
}

TEST_F(CombinedAirWindSkyBCTester, All)
{
  setupBCs({1, 2, 3, 4, 5}, {6, 7, 8, 9, 10});
  checkFlux();
  checkFluxDeriv();
  checkFluxDerivTair();
  checkFlux_rev();
}