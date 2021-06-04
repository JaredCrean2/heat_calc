#include "gtest/gtest.h"

#include "mesh/mesh_generator.h"
#include "mesh/mesh.h"
#include <memory>

TEST(Default, All)
{
  EXPECT_TRUE(true);
}

std::shared_ptr<Mesh::MeshCG> makeStandardMesh()
{
  Mesh::MeshSpec meshspec;
  meshspec.xmin = 0;
  meshspec.xmax = 2;
  meshspec.ymin = 0;
  meshspec.ymax = 2;
  meshspec.zmin = 0;
  meshspec.zmax = 2;
  meshspec.nx = 5;
  meshspec.ny = 5;
  meshspec.nz = 5;

  auto generator = Mesh::make_mesh_generator(meshspec, &(Mesh::identity));
  auto m = generator.generate();

  // make volume groups
  Mesh::MeshEntityGroupSpec volume_group("volume0");
  volume_group.addModelEntity(Mesh::ModelEntitySpec(3, 0));
  std::vector<Mesh::MeshEntityGroupSpec> volume_groups{volume_group};

  // make surface groups
  std::vector<Mesh::MeshEntityGroupSpec> surface_groups;
  for (int i=0; i < 6; ++i)
  {
    surface_groups.emplace_back(std::string("surface") + std::to_string(i));
    surface_groups.back().addModelEntity(Mesh::ModelEntitySpec(2, i), Mesh::ModelEntitySpec(3, 0));
    surface_groups.back().setIsDirichlet(true);
  }

  std::vector<Mesh::MeshEntityGroupSpec> other_surfaces;
  auto mesh = std::make_shared<Mesh::MeshCG>(m, volume_groups, surface_groups,
                                       other_surfaces, 1, 1);

  return mesh;
}

TEST(Mesh, Entities)
{
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


