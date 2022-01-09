
#include "gtest/gtest.h"
#include "discretization/DirichletBC_defs.h"
#include "discretization/disc_vector.h"
#include "discretization/source_term.h"
#include "physics/heat/HeatEquation.h"
#include "mesh_helper.h"

namespace {
  class HeatMMSTester : public StandardDiscSetup,
                        public testing::Test
  {
    protected:
      using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

      HeatMMSTester()
      {
        setup();
      }

      template <typename Tex, typename Tsrc>
      void setSolution(Tex ex_sol, Tsrc src)
      {
        heat = std::make_shared<Heat::HeatEquation>(disc);
        u_vec = makeDiscVector(disc);
        res_vec = makeDiscVector(disc);


        for (int i=0; i < disc->getNumVolDiscs(); ++i)
        {
          heat->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
          heat->addVolumeGroupParams(Heat::VolumeGroupParams(1, 1, 1));
        }

        for (int i=0; i < disc->getNumSurfDiscs(); ++i)
          heat->addDirichletBC(makeDirichletBCMMS(disc->getSurfDisc(i), ex_sol));

        res_vec->set(0);
        auto f = [&](Real x, Real y, Real z)
                    { return ex_sol(x, y, z, 0); };
        u_vec->setFunc(f);
      }

      HeatPtr heat;
      DiscVectorPtr u_vec;
      DiscVectorPtr res_vec;
    };

Real ex_sol(Real x, Real y, Real z, Real t, int degree)
{
  return std::pow(x, degree); // + std::pow(y, degree) + std::pow(z, degree);
}

Real src_func_dir(Real x, int degree)
{
  if (degree == 0 || degree == 1)
    return 0;
  else
    return -(degree * (degree - 1) * std::pow(x, degree-2));
}

Real src_func(Real x, Real y, Real z, Real t, int degree)
{
  return src_func_dir(x, degree); // + src_func_dir(y, degree) + src_func_dir(z, degree);
}

}  // namespace


TEST_F(HeatMMSTester, Constant)
{
  auto ex_sol = [] (Real x, Real y, Real z, Real t)
                   { return 1; };

  auto src_func = [] (Real x, Real y, Real z, Real t)
                     { return 0; };
  
  //setup();
  setSolution(ex_sol, src_func);

  heat->computeRhs(u_vec, 0.0, res_vec);
  res_vec->syncArrayToVector();

  auto& vec = res_vec->getVector();
  for (int i=0; i < vec.shape()[0]; ++i)
  {
    EXPECT_LE(std::abs(vec[i]), 1e-13);
  }
}


TEST_F(HeatMMSTester, PolynomialExactness)
{
  for (int sol_degree=1; sol_degree <= 3; ++sol_degree)
  {
    std::cout << "testing sol degree " << sol_degree << std::endl;
    for (int degree=0; degree <= sol_degree; ++degree)
    {
      std::cout << "  testing polynomial degree " << degree << std::endl;
      auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return ex_sol(x, y, z, t, degree); };

      auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                            { return src_func(x, y, z, t, degree); };

      setup(2*sol_degree, sol_degree);
      setSolution(ex_sol_l, src_func_l);

      heat->computeRhs(u_vec, 0.0, res_vec);
      res_vec->syncArrayToVector();

      auto& vec = res_vec->getVector();
      for (int i=0; i < vec.shape()[0]; ++i)
      {
        EXPECT_LE(std::abs(vec[i]), 1e-13);
      }
    }
  }
}