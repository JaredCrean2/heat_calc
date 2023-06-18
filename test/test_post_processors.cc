#include "gtest/gtest.h"
#include "discretization/NeumannBC_defs.h"
#include "discretization/disc_vector.h"
#include "mesh/mesh_input.h"
#include "physics/post_processors.h"
#include "test_helper.h"
#include "mesh_helper.h"
#include "discretization/surface_discretization.h"
#include <array>

namespace {
class PostProcessorTester : public ::testing::Test,
                            public StandardDiscSetup

{
  protected:
    PostProcessorTester()
    {
      setup(3, 3, Mesh::getMeshSpec(0, 1, 0, 2, 0, 3, 3, 3, 3), {false, false, false, false, false, false});
    }
};

}

TEST_F(PostProcessorTester, SurfaceAverage)
{
  auto u = makeDiscVector(disc);
  AuxiliaryEquationsStoragePtr u_aux = nullptr;
  u->set(2);

  auto surf = disc->getSurfDisc(1);

  auto f = [](double val) { return 2*val; };

  auto postproc = physics::makePostProcessorSurfaceIntegralAverage(surf, "postprocessor", f, MPI_COMM_WORLD);

  EXPECT_EQ(postproc->numValues(), 1);
  EXPECT_EQ(postproc->getNames()[0], "postprocessor");
  EXPECT_NEAR(postproc->getValues(u, u_aux, 0)[0], 4, 1e-13);
}

TEST_F(PostProcessorTester, BCFlux)
{
  auto u = makeDiscVector(disc);
  AuxiliaryEquationsStoragePtr u_aux = nullptr;
  u->set(2);

  auto surf = disc->getSurfDisc(0);

  auto f = [](double x, double y, double z, double t) { return std::array<Real, 3>{0, -2, 0}; };
  auto bc = makeNeumannBCMMS(surf, f);

  auto postproc = std::make_shared<physics::PostProcessorBCFlux>("postprocessor", bc, MPI_COMM_WORLD);


  EXPECT_EQ(postproc->numValues(), 1);
  EXPECT_EQ(postproc->getNames()[0], "postprocessor");
  EXPECT_NEAR(postproc->getValues(u, u_aux, 0)[0], 6, 1e-13);
}