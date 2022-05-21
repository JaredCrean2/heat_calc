#include "gtest/gtest.h"
#include "test_helper.h"
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
      SERIAL_ONLY_RETURN();

      m_meshspec1 = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 3, 4, 5);
      m_meshspec2 = Mesh::getMeshSpec(0, 1, 1, 2, 0, 1, 3, 5, 5);
      std::vector<Mesh::MeshSpec> meshspecs{m_meshspec1, m_meshspec2};

      auto generator = Mesh::make_mesh_generator(meshspecs, &Mesh::identity);
      m_mesh = generator.generate();
    }

    ~MeshGenerator2BlockTester() 
    {
      if (m_mesh)
        apf::destroyMesh(m_mesh);
    }

    apf::Mesh2* m_mesh = nullptr;
    Mesh::MeshSpec m_meshspec1;
    Mesh::MeshSpec m_meshspec2;
};

void countClassification(apf::Mesh2* mesh, ArrayType<int, 3>& counts)
{
  for (int dim=0; dim <= 3; ++dim)
    for (int tag=0; tag < counts.shape()[1]; ++tag)
      for (int dim2=0; dim <= 3; ++dim)
        counts[dim][tag][dim2] = 0;

  for (int dim=0; dim <= 3; ++dim)
  {
    apf::MeshIterator* it = mesh->begin(dim);
    apf::MeshEntity* e;
    while ( (e = mesh->iterate(it)) )
    {
      apf::ModelEntity* me = mesh->toModel(e);
      int me_dim = mesh->getModelType(me);
      int me_tag = mesh->getModelTag(me);
      counts[me_dim][me_tag][dim]++;
    }

    mesh->end(it);
  }
}


}

TEST_F(MeshGenerator2BlockTester, EntityCounts)
{
  SERIAL_ONLY();

  int nverts_x = m_meshspec1.nx + 1;
  int nverts_y = m_meshspec1.ny + m_meshspec2.ny + 1;
  int nverts_z = m_meshspec1.nz + 1;

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

TEST_F(MeshGenerator2BlockTester, EntityClassificationCounts)
{
  SERIAL_ONLY();
  
  ArrayType<int, 3> expected_counts(boost::extents[4][20][4]);  // geometric dimension x geometric tag x mesh entity dimension
  for (int i=0; i < 4; ++i)
    for (int j=0; j < 20; ++j)
      for (int k=0; k < 4; ++k)
        expected_counts[i][j][k] = 0;

  for (int i=0; i < 12; ++i)
    expected_counts[0][i][0] = 1;

  // mesh verts classified on geometric edges
  expected_counts[1][0][0]  = m_meshspec1.nx - 1;
  expected_counts[1][1][0]  = m_meshspec1.ny - 1;
  expected_counts[1][2][0]  = m_meshspec1.nx - 1;
  expected_counts[1][3][0]  = m_meshspec1.ny - 1;
  expected_counts[1][4][0]  = m_meshspec1.nx - 1;
  expected_counts[1][5][0]  = m_meshspec1.ny - 1;
  expected_counts[1][6][0]  = m_meshspec1.nx - 1;
  expected_counts[1][7][0]  = m_meshspec1.ny - 1;
  expected_counts[1][8][0]  = m_meshspec1.nz - 1;
  expected_counts[1][9][0]  = m_meshspec1.nz - 1;
  expected_counts[1][10][0] = m_meshspec1.nz - 1;
  expected_counts[1][11][0] = m_meshspec1.nz - 1;

  expected_counts[1][12][0] = m_meshspec2.ny - 1;
  expected_counts[1][13][0] = m_meshspec2.nx - 1;
  expected_counts[1][14][0] = m_meshspec2.ny - 1;
  expected_counts[1][15][0] = m_meshspec2.ny - 1;
  expected_counts[1][16][0] = m_meshspec2.nx - 1;
  expected_counts[1][17][0] = m_meshspec2.ny - 1;
  expected_counts[1][18][0] = m_meshspec2.nz - 1;
  expected_counts[1][19][0] = m_meshspec2.nz - 1;

  // mesh edges classified on geometric edges
  expected_counts[1][0][1]  = m_meshspec1.nx ;
  expected_counts[1][1][1]  = m_meshspec1.ny ;
  expected_counts[1][2][1]  = m_meshspec1.nx ;
  expected_counts[1][3][1]  = m_meshspec1.ny ;
  expected_counts[1][4][1]  = m_meshspec1.nx ;
  expected_counts[1][5][1]  = m_meshspec1.ny ;
  expected_counts[1][6][1]  = m_meshspec1.nx ;
  expected_counts[1][7][1]  = m_meshspec1.ny ;
  expected_counts[1][8][1]  = m_meshspec1.nz ;
  expected_counts[1][9][1]  = m_meshspec1.nz ;
  expected_counts[1][10][1] = m_meshspec1.nz ;
  expected_counts[1][11][1] = m_meshspec1.nz ;

  expected_counts[1][12][1] = m_meshspec2.ny;
  expected_counts[1][13][1] = m_meshspec2.nx;
  expected_counts[1][14][1] = m_meshspec2.ny;
  expected_counts[1][15][1] = m_meshspec2.ny;
  expected_counts[1][16][1] = m_meshspec2.nx;
  expected_counts[1][17][1] = m_meshspec2.ny;
  expected_counts[1][18][1] = m_meshspec2.nz;
  expected_counts[1][19][1] = m_meshspec2.nz;

  // mesh verts classified on model faces
  expected_counts[2][0][0]  = (m_meshspec1.nx - 1)*(m_meshspec1.nz - 1);
  expected_counts[2][1][0]  = (m_meshspec1.ny - 1)*(m_meshspec1.nz - 1);
  expected_counts[2][2][0]  = (m_meshspec1.nx - 1)*(m_meshspec1.nz - 1);
  expected_counts[2][3][0]  = (m_meshspec1.ny - 1)*(m_meshspec1.nz - 1);
  expected_counts[2][4][0]  = (m_meshspec1.nx - 1)*(m_meshspec1.ny - 1);
  expected_counts[2][5][0]  = (m_meshspec1.nx - 1)*(m_meshspec1.ny - 1);

  expected_counts[2][6][0]  = (m_meshspec2.ny - 1)*(m_meshspec2.nz - 1);
  expected_counts[2][7][0]  = (m_meshspec2.nx - 1)*(m_meshspec2.nz - 1);
  expected_counts[2][8][0]  = (m_meshspec2.ny - 1)*(m_meshspec2.nz - 1);
  expected_counts[2][9][0]  = (m_meshspec2.nx - 1)*(m_meshspec2.ny - 1);
  expected_counts[2][10][0] = (m_meshspec2.nx - 1)*(m_meshspec2.ny - 1);
  
  // mesh edges classified on model faces
  expected_counts[2][0][1]  = (m_meshspec1.nx - 1) * m_meshspec1.nz + (m_meshspec1.nz - 1) * m_meshspec1.nx;
  expected_counts[2][1][1]  = (m_meshspec1.ny - 1) * m_meshspec1.nz + (m_meshspec1.nz - 1) * m_meshspec1.ny;
  expected_counts[2][2][1]  = (m_meshspec1.nx - 1) * m_meshspec1.nz + (m_meshspec1.nz - 1) * m_meshspec1.nx;
  expected_counts[2][3][1]  = (m_meshspec1.ny - 1) * m_meshspec1.nz + (m_meshspec1.nz - 1) * m_meshspec1.ny;
  expected_counts[2][4][1]  = (m_meshspec1.nx - 1) * m_meshspec1.ny + (m_meshspec1.ny - 1) * m_meshspec1.nx;
  expected_counts[2][5][1]  = (m_meshspec1.nx - 1) * m_meshspec1.ny + (m_meshspec1.ny - 1) * m_meshspec1.nx;

  expected_counts[2][6][1]  = (m_meshspec2.ny - 1) * m_meshspec2.nz + (m_meshspec2.nz - 1) * m_meshspec2.ny;
  expected_counts[2][7][1]  = (m_meshspec2.nx - 1) * m_meshspec2.nz + (m_meshspec2.nz - 1) * m_meshspec2.nx;
  expected_counts[2][8][1]  = (m_meshspec2.ny - 1) * m_meshspec2.nz + (m_meshspec2.nz - 1) * m_meshspec2.ny;
  expected_counts[2][9][1]  = (m_meshspec2.nx - 1) * m_meshspec2.nz + (m_meshspec2.nz - 1) * m_meshspec2.nx;
  expected_counts[2][10][1] = (m_meshspec2.nx - 1) * m_meshspec2.nz + (m_meshspec2.nz - 1) * m_meshspec2.nx;

  // mesh faces classified on model faces
  expected_counts[2][0][2]  = m_meshspec1.nx * m_meshspec1.nz;
  expected_counts[2][1][2]  = m_meshspec1.ny * m_meshspec1.nz;
  expected_counts[2][2][2]  = m_meshspec1.nx * m_meshspec1.nz;
  expected_counts[2][3][2]  = m_meshspec1.ny * m_meshspec1.nz;
  expected_counts[2][4][2]  = m_meshspec1.nx * m_meshspec1.ny;
  expected_counts[2][5][2]  = m_meshspec1.nx * m_meshspec1.ny;

  expected_counts[2][6][2]  = m_meshspec2.ny * m_meshspec2.nz;
  expected_counts[2][7][2]  = m_meshspec2.nx * m_meshspec2.nz;
  expected_counts[2][8][2]  = m_meshspec2.ny * m_meshspec2.nz;
  expected_counts[2][9][2]  = m_meshspec2.nx * m_meshspec2.ny;
  expected_counts[2][10][2] = m_meshspec2.nx * m_meshspec2.ny;


  // mesh verts classified on model regions
  expected_counts[3][0][0]  = (m_meshspec1.nx - 1) * (m_meshspec1.ny - 1) * (m_meshspec1.nz - 1);
  expected_counts[3][1][0]  = (m_meshspec2.nx - 1) * (m_meshspec2.ny - 1) * (m_meshspec2.nz - 1);

  // mesh edges classified on model regions
  expected_counts[3][0][1] = (m_meshspec1.nx - 1) * (m_meshspec1.ny - 1) *  m_meshspec1.nz      + 
                              m_meshspec1.nx      * (m_meshspec1.ny - 1) * (m_meshspec1.nz - 1) +
                             (m_meshspec1.nx - 1) *  m_meshspec1.ny      * (m_meshspec1.nz - 1);
  expected_counts[3][1][1] = (m_meshspec2.nx - 1) * (m_meshspec2.ny - 1) *  m_meshspec2.nz      + 
                              m_meshspec2.nx      * (m_meshspec2.ny - 1) * (m_meshspec2.nz - 1) +
                             (m_meshspec2.nx - 1) *  m_meshspec2.ny      * (m_meshspec2.nz - 1);

  // mesh faces classified on model regions
  expected_counts[3][0][2] = m_meshspec1.nx * m_meshspec1.ny * (m_meshspec1.nz - 1) +
                             m_meshspec1.ny * m_meshspec1.nz * (m_meshspec1.nx - 1) +
                             m_meshspec1.nz * m_meshspec1.nx * (m_meshspec1.ny - 1);
  expected_counts[3][1][2] = m_meshspec2.nx * m_meshspec2.ny * (m_meshspec2.nz - 1) +
                             m_meshspec2.ny * m_meshspec2.nz * (m_meshspec2.nx - 1) +
                             m_meshspec2.nz * m_meshspec2.nx * (m_meshspec2.ny - 1);

  // mesh elements classified on model elements
  expected_counts[3][0][3] = m_meshspec1.nx * m_meshspec1.ny * m_meshspec1.nz;
  expected_counts[3][1][3] = m_meshspec2.nx * m_meshspec2.ny * m_meshspec2.nz;

  ArrayType<int, 3> actual_counts(boost::extents[4][20][4]);
  countClassification(m_mesh, actual_counts);

  for (int me_dim=0; me_dim <= 3; ++me_dim)
    for (int tag=0; tag < 20; ++tag)
      for (int mesh_dim=0; mesh_dim <= 3; ++mesh_dim)
      {
        EXPECT_EQ(actual_counts[me_dim][tag][mesh_dim], expected_counts[me_dim][tag][mesh_dim]);
      }


}