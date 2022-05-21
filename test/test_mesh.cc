#include "gtest/gtest.h"
#include "test_helper.h"
#include "mesh_helper.h"
#include "mesh/mesh.h"

TEST(Default, All)
{
  EXPECT_TRUE(true);
}

TEST(Mesh, Entities)
{
  SERIAL_ONLY();

  auto mesh = makeStandardMesh();
  EXPECT_EQ(mesh->getNumBCSurfaces(), 6);
  EXPECT_EQ(mesh->getNumSurfaces(), 6);
  EXPECT_EQ(mesh->getNumVolumeGroups(), 1);
  EXPECT_EQ(mesh->getNumTotalDofs(), 6*6*6);
  // number of dofs in cube - dofs on interior of faces
  // - dofs on interior of edge - dofs on vertices
  int ndofs = 6*6*6 - 4*4*6 - (6-2)*12 - 8;
  EXPECT_EQ(mesh->getNumDofs(), ndofs);

  auto& vol_group = mesh->getElements(0);
  EXPECT_EQ(vol_group.nodenums.shape()[0], std::size_t(5*5*5));

  for (std::size_t i=0; i < vol_group.nodenums.shape()[0]; ++i)
    for (std::size_t j=0; j < vol_group.nodenums.shape()[1]; ++j)
    {
      EXPECT_TRUE(vol_group.nodenums[i][j] >= 0);
      EXPECT_TRUE(vol_group.nodenums[i][j] < mesh->getNumTotalDofs());
    }

  for (int i=0; i < 6; ++i)
  {
    auto& face_group = mesh->getFaces(i);
    EXPECT_EQ(face_group.faces.size(), std::size_t(5*5));
  }
}


