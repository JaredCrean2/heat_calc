
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
}


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
    std::cout << "vec " << i << " = " << vec[i] << std::endl;
    EXPECT_LE(std::abs(vec[i]), 1e-13);
  }
}