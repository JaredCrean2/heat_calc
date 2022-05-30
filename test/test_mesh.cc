#include "gtest/gtest.h"
#include <PCU.h>
#include "test_helper.h"
#include "mesh_helper.h"
#include "mesh/mesh.h"

namespace {
  void test_all_dofs_present(std::shared_ptr<Mesh::MeshCG> mesh)
  {
    int num_owned_dofs = mesh->getNumOwnedDofs(), num_global_dofs;
    MPI_Allreduce(&num_owned_dofs, &num_global_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    std::vector<int> used_global_dofs(num_global_dofs, false);
    std::vector<DofInt> owned_to_local_dofs, local_to_global_dofs;
    mesh->getOwnedLocalDofInfo(owned_to_local_dofs);
    mesh->getLocalToGlobalDofs(local_to_global_dofs);

    for (int i=0; i < num_owned_dofs; ++i)
    {
      auto local_dof = owned_to_local_dofs[i];
      auto global_dof = local_to_global_dofs[local_dof];
      used_global_dofs[global_dof] = true;
    }

    std::vector<int> used_global_dofs_output(num_global_dofs, false);
    MPI_Allreduce(used_global_dofs.data(), used_global_dofs_output.data(), num_global_dofs, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    for (int i=0; i < num_global_dofs; ++i)
      EXPECT_TRUE(used_global_dofs_output[i]);
  }
}


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


TEST(Mesh, ParallelEntities)
{
  auto mesh = makeStandardMesh();

  EXPECT_EQ(mesh->getNumBCSurfaces(), 6);
  EXPECT_EQ(mesh->getNumSurfaces(), 6);
  EXPECT_EQ(mesh->getNumVolumeGroups(), 1);

  int num_owned_dofs = mesh->getNumOwnedDofs(), num_global_dofs;
  MPI_Allreduce(&num_owned_dofs, &num_global_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  EXPECT_EQ(num_global_dofs, 4*4*4);

  test_all_dofs_present(mesh);
/*
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
*/
}

