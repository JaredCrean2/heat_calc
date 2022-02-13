#include "gtest/gtest.h"
#include "mesh/reference_element_def.h"
#include "mesh/reference_element_apf.h"


TEST(ReferenceElementApf, Linear)
{
  auto ref_el = reference_element::getLagrangeHexReferenceElement(1);
  auto fshape = apf::getHexFieldShape(ref_el);

  EXPECT_TRUE(fshape->hasNodesIn(0));
  EXPECT_FALSE(fshape->hasNodesIn(1));
  EXPECT_FALSE(fshape->hasNodesIn(2));
  EXPECT_FALSE(fshape->hasNodesIn(3));

  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::VERTEX), 1);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::EDGE),   0);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::QUAD),   0);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::HEX),    0);

  auto vertshape = fshape->getEntityShape(apf::Mesh::VERTEX);
  EXPECT_EQ(vertshape->countNodes(), 1);

  auto edgeshape = fshape->getEntityShape(apf::Mesh::EDGE);
  EXPECT_EQ(edgeshape->countNodes(), 2);

  auto quadshape = fshape->getEntityShape(apf::Mesh::QUAD);
  EXPECT_EQ(quadshape->countNodes(), 4);

  auto hexshape = fshape->getEntityShape(apf::Mesh::HEX);
  EXPECT_EQ(hexshape->countNodes(), 8);
}


TEST(ReferenceElementApf, Quadratic)
{
  auto ref_el = reference_element::getLagrangeHexReferenceElement(2);
  auto fshape = apf::getHexFieldShape(ref_el);

  EXPECT_TRUE(fshape->hasNodesIn(0));
  EXPECT_TRUE(fshape->hasNodesIn(1));
  EXPECT_TRUE(fshape->hasNodesIn(2));
  EXPECT_TRUE(fshape->hasNodesIn(3));

  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::VERTEX), 1);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::EDGE),   1);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::QUAD),   1);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::HEX),    1);

  auto vertshape = fshape->getEntityShape(apf::Mesh::VERTEX);
  EXPECT_EQ(vertshape->countNodes(), 1);

  auto edgeshape = fshape->getEntityShape(apf::Mesh::EDGE);
  EXPECT_EQ(edgeshape->countNodes(), 3);

  auto quadshape = fshape->getEntityShape(apf::Mesh::QUAD);
  EXPECT_EQ(quadshape->countNodes(), 9);

  auto hexshape = fshape->getEntityShape(apf::Mesh::HEX);
  EXPECT_EQ(hexshape->countNodes(), 27);
}


TEST(ReferenceElementApf, Cubic)
{
  auto ref_el = reference_element::getLagrangeHexReferenceElement(3);
  auto fshape = apf::getHexFieldShape(ref_el);

  EXPECT_TRUE(fshape->hasNodesIn(0));
  EXPECT_TRUE(fshape->hasNodesIn(1));
  EXPECT_TRUE(fshape->hasNodesIn(2));
  EXPECT_TRUE(fshape->hasNodesIn(3));

  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::VERTEX), 1);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::EDGE),   2);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::QUAD),   4);
  EXPECT_EQ(fshape->countNodesOn(apf::Mesh::HEX),    8);

  auto vertshape = fshape->getEntityShape(apf::Mesh::VERTEX);
  EXPECT_EQ(vertshape->countNodes(), 1);

  auto edgeshape = fshape->getEntityShape(apf::Mesh::EDGE);
  EXPECT_EQ(edgeshape->countNodes(), 4);

  auto quadshape = fshape->getEntityShape(apf::Mesh::QUAD);
  EXPECT_EQ(quadshape->countNodes(), 16);

  auto hexshape = fshape->getEntityShape(apf::Mesh::HEX);
  EXPECT_EQ(hexshape->countNodes(), 64);
}