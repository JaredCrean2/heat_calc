#include "gtest/gtest.h"
#include "mesh/mesh_geometry_multi.h"
#include "mesh/mesh_generator_multi_block.h"

namespace {

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

double computeMeshVolume(apf::Mesh* m)
{
  double volume = 0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(3);
  while ( (e = m->iterate(it)) )
    volume += apf::measure(m, e);
  m->end(it);

  return volume;
}

std::map<int, double> computeMeshVolumePerDomain(apf::Mesh* m)
{
  std::map<int, double> volumes;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(3);
  while ( (e = m->iterate(it)) )
  {
    apf::ModelEntity* me = m->toModel(e);
    int tag = m->getModelTag(me);

    volumes[tag] += apf::measure(m, e);
  }
  m->end(it);  

  return volumes;
}

}

TEST(MeshGeneratorGeneral, MiddleBlock)
{
  Mesh::MultiBlockMeshSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2, 0, 3, 0, 4, 3, 4, 5);
  spec.create_middle_block = true;

  Mesh::MeshGeneratorMultiBlock generator(spec);
  apf::Mesh2* m = generator.generate();
  EXPECT_NEAR(computeMeshVolume(m), 2*3*4, 1e-13);
}

TEST(MeshGeneratorGeneral, MiddleBlockPlusOneLayerY)
{
  Mesh::MultiBlockMeshSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2, 0, 3, 0, 4, 3, 4, 5);
  spec.create_middle_block = true;
  spec.numel_plusy.push_back(6);
  spec.thickness_plusy.push_back(2);

  Mesh::MeshGeneratorMultiBlock generator(spec);
  apf::Mesh2* m = generator.generate();
  EXPECT_NEAR(computeMeshVolume(m), 2*3*4 + 2*4*2, 1e-13);
}

TEST(MeshGeneratorGeneral, MiddleBlockPlusTwoLayerY)
{
  Mesh::MultiBlockMeshSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2, 0, 3, 0, 4, 3, 4, 5);
  spec.create_middle_block = true;
  spec.numel_plusy.push_back(6);     spec.numel_plusy.push_back(7);
  spec.thickness_plusy.push_back(2); spec.thickness_plusy.push_back(3);

  Mesh::MeshGeneratorMultiBlock generator(spec);
  apf::Mesh2* m = generator.generate();
  EXPECT_NEAR(computeMeshVolume(m), 2*3*4 + 2*4*2 + 2*4*3, 1e-12);
}


TEST(MeshGeneratorGeneral, MiddleBlockPlusTwoLayers)
{
  Mesh::MultiBlockMeshSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2, 0, 3, 0, 4, 3, 4, 5);
  spec.create_middle_block = true;

  spec.numel_plusx.push_back(3);     spec.numel_plusx.push_back(3);
  spec.thickness_plusx.push_back(2); spec.thickness_plusx.push_back(3);

  spec.numel_plusy.push_back(3);     spec.numel_plusy.push_back(3);
  spec.thickness_plusy.push_back(2); spec.thickness_plusy.push_back(3);

  spec.numel_plusz.push_back(3);     spec.numel_plusz.push_back(3);
  spec.thickness_plusz.push_back(2); spec.thickness_plusz.push_back(3);

  Mesh::MeshGeneratorMultiBlock generator(spec);
  apf::Mesh2* m = generator.generate();
  EXPECT_NEAR(computeMeshVolume(m), (2+2+3)*(3+2+3)*(4+2+3), 1e-11);

  auto volumes = computeMeshVolumePerDomain(m);
  for (auto p : volumes)
    std::cout << "domain " << p.first << " has volume " << p.second << std::endl;
}