#include "gtest/gtest.h"
#include "test_helper.h"
#include "mesh/reference_element.h"

namespace 
{
bool checkTensorProductXiRange(Mesh::ReferenceElement* ref_el)
{
  auto& node_xi = ref_el->getTensorProductXi();
  auto range = ref_el->getXiRange();
  bool is_in_range = true;

  for (auto xi : node_xi)
    is_in_range = is_in_range && xi >= range.first && xi <= range.second;

  return is_in_range;
}

bool checkTensorProductXiUniqueness(Mesh::ReferenceElement* ref_el)
{
  auto& node_xi = ref_el->getTensorProductXi();
  bool is_unique = true;

  for (unsigned int i=0; i < node_xi.size(); ++i)
    for (unsigned int j=0; j < node_xi.size(); ++j)
      if (i != j)
        is_unique = is_unique && std::abs(node_xi[i] - node_xi[j]) > 1e-13;

  return is_unique;
}

bool checkTensorProductNodeMapRange(Mesh::ReferenceElement* ref_el)
{
  auto& nodemap = ref_el->getTensorProductMap();
  for (int d=0; d < 3; ++d)
    EXPECT_EQ(nodemap.shape()[d], ref_el->getNumNodesTensorProduct());

  bool is_in_range = true;
  for (int i=0; i < ref_el->getNumNodesTensorProduct(); ++i)
    for (int j=0; j < ref_el->getNumNodesTensorProduct(); ++j)
      for (int k=0; k < ref_el->getNumNodesTensorProduct(); ++k)
        is_in_range = is_in_range && nodemap[i][j][k] >= 0
                                  && nodemap[i][j][k] < ref_el->getNumNodes();
  
  return is_in_range;
}

bool checkTensorProductNodeMapUniqueness(Mesh::ReferenceElement* ref_el)
{
  auto& nodemap = ref_el->getTensorProductMap();
  std::set<LocalIndex> nodes;
  for (int i=0; i < ref_el->getNumNodesTensorProduct(); ++i)
    for (int j=0; j < ref_el->getNumNodesTensorProduct(); ++j)
      for (int k=0; k < ref_el->getNumNodesTensorProduct(); ++k)
        nodes.insert(nodemap[i][j][k]);

  
  return nodes.size() == ref_el->getNumNodes();
}

}

TEST(ReferenceElement, Counts)
{
  SERIAL_ONLY();

  Mesh::ReferenceElement* hex1 = Mesh::getReferenceElement(apf::Mesh::HEX, 1);
  Mesh::ReferenceElement* hex2 = Mesh::getReferenceElement(apf::Mesh::HEX, 2);
  Mesh::ReferenceElement* hex3 = Mesh::getReferenceElement(apf::Mesh::HEX, 3);

  EXPECT_EQ(hex1->getDegree(), 1);
  EXPECT_EQ(hex2->getDegree(), 2);
  EXPECT_EQ(hex3->getDegree(), 3);

  EXPECT_EQ(hex1->getNumFaces(), 6);
  EXPECT_EQ(hex2->getNumFaces(), 6);
  EXPECT_EQ(hex3->getNumFaces(), 6);

  EXPECT_EQ(hex1->getNumNodes(), 8);
  EXPECT_EQ(hex2->getNumNodes(), 27);
  EXPECT_EQ(hex3->getNumNodes(), 64);

  EXPECT_EQ(hex1->getNumNodesTensorProduct(), 2);
  EXPECT_EQ(hex2->getNumNodesTensorProduct(), 3);
  EXPECT_EQ(hex3->getNumNodesTensorProduct(), 4);
}


TEST(ReferenceElement, getTensorProductXi)
{
  SERIAL_ONLY();

  Mesh::ReferenceElement* hex1 = Mesh::getReferenceElement(apf::Mesh::HEX, 1);
  Mesh::ReferenceElement* hex2 = Mesh::getReferenceElement(apf::Mesh::HEX, 2);
  Mesh::ReferenceElement* hex3 = Mesh::getReferenceElement(apf::Mesh::HEX, 3);

  EXPECT_TRUE(checkTensorProductXiRange(hex1));
  EXPECT_TRUE(checkTensorProductXiRange(hex2));
  EXPECT_TRUE(checkTensorProductXiRange(hex3));

  EXPECT_TRUE(checkTensorProductXiUniqueness(hex1));
  EXPECT_TRUE(checkTensorProductXiUniqueness(hex2));
  EXPECT_TRUE(checkTensorProductXiUniqueness(hex3));
}

TEST(ReferenceElement, getTensorProductMap)
{
  SERIAL_ONLY();
  
  Mesh::ReferenceElement* hex1 = Mesh::getReferenceElement(apf::Mesh::HEX, 1);
  Mesh::ReferenceElement* hex2 = Mesh::getReferenceElement(apf::Mesh::HEX, 2);
  Mesh::ReferenceElement* hex3 = Mesh::getReferenceElement(apf::Mesh::HEX, 3);

  EXPECT_TRUE(checkTensorProductNodeMapRange(hex1));
  EXPECT_TRUE(checkTensorProductNodeMapRange(hex2));
  EXPECT_TRUE(checkTensorProductNodeMapRange(hex3));


  EXPECT_TRUE(checkTensorProductNodeMapUniqueness(hex1));
  EXPECT_TRUE(checkTensorProductNodeMapUniqueness(hex2));
  EXPECT_TRUE(checkTensorProductNodeMapUniqueness(hex3));
}