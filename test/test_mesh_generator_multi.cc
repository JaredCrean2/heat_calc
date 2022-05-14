#include "gtest/gtest.h"
#include <apf.h>

#include "mesh/mesh_generator.h"
#include "mesh/mesh_generator_multi.h"
#include "mesh/mesh_input.h"

namespace {

class MeshGenerator2BlockTester : public ::testing::Test
{
  protected:
    MeshGenerator2BlockTester()
    {
      /*
      m_meshspec1 = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 3, 4, 5);
      m_meshspec2 = m_meshspec1;
      m_meshspec2.ymin = m_meshspec1.ymax;
      m_meshspec2.ymax = 2;
      m_meshspec2.ny = 6;

      auto generator = Mesh::make_mesh_generator({m_meshspec1, m_meshspec2}, &Mesh::identity);
      m_mesh = generator.generate();
      */

      m_meshspec1 = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 3, 4, 5);
      m_meshspec2 = Mesh::getMeshSpec(0, 0, 0, 0, 0, 0, 0, 0, 0);
      std::vector<Mesh::MeshSpec> meshspecs{m_meshspec1};

      auto generator = Mesh::make_mesh_generator(meshspecs, &Mesh::identity);
      m_mesh = generator.generate();
    }

    ~MeshGenerator2BlockTester() {apf::destroyMesh(m_mesh); }

    apf::Mesh2* m_mesh;
    Mesh::MeshSpec m_meshspec1;
    Mesh::MeshSpec m_meshspec2;
};


}

TEST_F(MeshGenerator2BlockTester, EntityCounts)
{
  int nverts_x = m_meshspec1.nx + m_meshspec2.nx + 1;
  int nverts_y = m_meshspec1.ny + m_meshspec2.ny + 1;
  int nverts_z = m_meshspec1.nz + m_meshspec2.nz + 1;

  int nedges_x = nverts_x - 1;
  int nedges_y = nverts_y - 1;
  int nedges_z = nverts_z - 1;
  int nedges_total = nverts_x * nverts_y * nedges_z +
                     nverts_x * nedges_y * nverts_z +
                     nedges_x * nverts_y * nverts_z;

  int nfaces_total = nedges_x * nedges_y * nverts_z +
                     nedges_x * nverts_y * nedges_z +
                     nverts_x * nedges_y * nedges_z;

  int nelements_total = nedges_x * nedges_y * nedges_z;


  EXPECT_EQ(m_mesh->count(0), nverts_x * nverts_y * nverts_z);
  EXPECT_EQ(m_mesh->count(1), nedges_total);
  EXPECT_EQ(m_mesh->count(2), nfaces_total);
  EXPECT_EQ(m_mesh->count(3), nelements_total);
}