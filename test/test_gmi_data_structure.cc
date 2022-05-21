#include "gtest/gtest.h"
#include "test_helper.h"
#include "mesh/gmiDataStructure.h"

namespace {
  void test_entity(const mesh_gmi::GMIEntity& entity, int dim, int tag, int index)
  {
    EXPECT_EQ(entity.getDim(), dim);
    EXPECT_EQ(entity.getTag(), tag);
    EXPECT_EQ(entity.getIndex(), index);
  }


  void test_set_and_free(struct gmi_model* gmi_model, struct gmi_set* adjacent, int expected_dim, const std::vector<int>& expected_tags)
  {
    EXPECT_EQ(adjacent->n, expected_tags.size());
    for (unsigned int i=0; i < expected_tags.size(); ++i)
    {
      EXPECT_EQ(gmi_dim(gmi_model, adjacent->e[i]), expected_dim);;
      EXPECT_EQ(gmi_tag(gmi_model, adjacent->e[i]), expected_tags[i]);
    }
    gmi_free_set(adjacent);
  }
}

TEST(GMIDataStructure, EntityRetrieval)
{
  SERIAL_ONLY();

  mesh_gmi::GMITopo topo;

  for (int dim=0; dim <= 3; ++dim)
  {
    topo.createEntity(dim, 2);
    topo.createEntity(dim, 3);
    topo.createEntity(dim, 4);
    EXPECT_EQ(topo.getNumEntities(dim), 3);

    test_entity(topo.getEntityByTag(dim, 2), dim, 2, 0);
    test_entity(topo.getEntityByTag(dim, 3), dim, 3, 1);
    test_entity(topo.getEntityByTag(dim, 4), dim, 4, 2);

    test_entity(topo.getEntityByIndex(dim, 0), dim, 2, 0);
    test_entity(topo.getEntityByIndex(dim, 1), dim, 3, 1);
    test_entity(topo.getEntityByIndex(dim, 2), dim, 4, 2);

    test_entity(*(topo.gmiFind(dim, 2)), dim, 2, 0);
    test_entity(*(topo.gmiFind(dim, 3)), dim, 3, 1);
    test_entity(*(topo.gmiFind(dim, 4)), dim, 4, 2);

    EXPECT_ANY_THROW(topo.getEntityByTag(dim, 0));
  }
}

TEST(GMIDataStructure, Iteration)
{
  SERIAL_ONLY();

  mesh_gmi::GMITopo topo;
  int dim=1;
  topo.createEntity(dim, 2);
  topo.createEntity(dim, 3);
  topo.createEntity(dim, 4);

  auto it = topo.gmiBegin(dim);
  mesh_gmi::GMIEntity* entity;
  int idx=0;
  while ( (entity = topo.gmiNext(it)))
  {
    test_entity(*entity, dim, idx+2, idx);
    idx++;
  }
  topo.gmiEnd(it);
}

TEST(GMIDataStructure, Adjacency)
{
  SERIAL_ONLY();
  
  auto topo = std::make_shared<mesh_gmi::GMITopo>();

  {
    auto& v1 = topo->createEntity(0, 1);
    auto& v2 = topo->createEntity(0, 2);

    auto& e1 = topo->createEntity(1, 1);
    auto& e2 = topo->createEntity(1, 2);
    auto& e3 = topo->createEntity(1, 3);

    auto& f1 = topo->createEntity(2, 1);
    auto& f2 = topo->createEntity(2, 2);

    v1.addUpward(e1);
    v2.addUpward(e1);
    e1.addDownward(v1);
    e1.addDownward(v2);

    v2.addUpward(e3);
    e3.addDownward(v2);

    e1.addUpward(f1);
    e2.addUpward(f1);
    e3.addUpward(f1);
    f1.addDownward(e1);
    f1.addDownward(e2);
    f1.addDownward(e3);

    e2.addUpward(f2);
    e3.addUpward(f2);
    f2.addDownward(e2);
    f2.addDownward(e3);
  }


  auto gmi_model = mesh_gmi::createGMITopo(topo);
  auto v1 = gmi_find(gmi_model, 0, 1);
  auto v2 = gmi_find(gmi_model, 0, 2);
  
  auto e1 = gmi_find(gmi_model, 1, 1);
  auto e2 = gmi_find(gmi_model, 1, 2);
  auto e3 = gmi_find(gmi_model, 1, 3);

  
  auto f1 = gmi_find(gmi_model, 2, 1);
  auto f2 = gmi_find(gmi_model, 2, 2);
  

  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, v1, 1), 1, {1});
  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, v2, 1), 1, {1, 3});

  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, e1, 0), 0, {1, 2});
  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, e1, 2), 2, {1});

  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, e2, 0), 0, {});
  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, e2, 2), 2, {1, 2});

  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, e3, 0), 0, {2});
  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, e3, 2), 2, {1, 2});

  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, f1, 1), 1, {1, 2, 3});
  test_set_and_free(gmi_model, gmi_adjacent(gmi_model, f2, 1), 1, {2, 3});

  EXPECT_NO_THROW(mesh_gmi::verify(*topo));

  gmi_model->ops->destroy(gmi_model);  
}