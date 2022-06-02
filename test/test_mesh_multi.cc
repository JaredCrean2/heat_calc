#include "mesh_helper.h"
#include "gtest/gtest.h"
#include "test_helper.h"

namespace {
  class MeshMultiTester : public ::testing::Test, public StandardDiscSetupMulti
  {};
}

TEST_F(MeshMultiTester, Counts)
{
  SERIAL_ONLY();

  setup(3, 1);

  EXPECT_EQ(mesh->getNumBCSurfaces(), 10);
  EXPECT_EQ(mesh->getNumSurfaces(), 11);
  EXPECT_EQ(mesh->getNumVolumeGroups(), 2);
}


TEST_F(MeshMultiTester, DofCounts)
{
  SERIAL_ONLY();
  
  setup(3, 1);

  int nverts_x = specs[0].nx + 1;
  int nverts_y = specs[0].ny + specs[1].ny + 1;
  int nverts_z = specs[0].nz + 1;

  int num_active_dofs = (nverts_x - 2) * (nverts_y - 2) * (nverts_z - 2);
  int num_dirichlet_dofs = nverts_x * nverts_y * nverts_z - num_active_dofs;

  EXPECT_EQ(num_active_dofs, mesh->getNumDofs());
  EXPECT_EQ(nverts_x*nverts_y*nverts_z, mesh->getNumTotalDofs());
  EXPECT_EQ(num_dirichlet_dofs, mesh->getNumTotalDofs() - mesh->getNumDofs());
  EXPECT_EQ(disc->getDofNumbering()->getNumDofs(), num_active_dofs);
}