#include "gtest/gtest.h"
#include "mesh/reference_element_interface.h"
#include "mesh/reference_element_geometry_interface.h"


TEST(ReferenceElementHex, Geometry)
{
  auto ref_el = reference_element::getLagrangeHexReferenceElement(1);

  EXPECT_EQ(ref_el->getNumVerts(), 8);
  EXPECT_EQ(ref_el->getNumEntities(0), 8);
  EXPECT_EQ(ref_el->getNumEdges(), 12);
  EXPECT_EQ(ref_el->getNumEntities(1), 12);
  EXPECT_EQ(ref_el->getNumFaces(), 6);
  EXPECT_EQ(ref_el->getNumEntities(2), 6);

  for (int i=0; i < ref_el->getNumVerts(); ++i)
    EXPECT_EQ(ref_el->getVert(i)->getDimension(), 0);

  for (int i=0; i < ref_el->getNumEdges(); ++i)
    EXPECT_EQ(ref_el->getEdge(i)->getDimension(), 1);

  for (int i=0; i < ref_el->getNumFaces(); ++i)
    EXPECT_EQ(ref_el->getFace(i)->getDimension(), 2);
}


TEST(ReferenceElementHex, Linear)
{
  auto ref_el = reference_element::getLagrangeHexReferenceElement(1);

  EXPECT_EQ(ref_el->getNumNodes(0), 1);
  EXPECT_EQ(ref_el->getNumNodes(1), 0);
  EXPECT_EQ(ref_el->getNumNodes(2), 0);
  EXPECT_EQ(ref_el->getNumNodes(3), 0);
  EXPECT_EQ(ref_el->getNumNodesTotal(), 8);

  for (int i=0; i < 8; ++i)
    EXPECT_EQ(ref_el->getNodeIndex(0, i, 0), i);
}


TEST(ReferenceElementHex, Quadratic)
{
  auto ref_el = reference_element::getLagrangeHexReferenceElement(2);

  EXPECT_EQ(ref_el->getNumNodes(0), 1);
  EXPECT_EQ(ref_el->getNumNodes(1), 1);
  EXPECT_EQ(ref_el->getNumNodes(2), 1);
  EXPECT_EQ(ref_el->getNumNodes(3), 1);
  EXPECT_EQ(ref_el->getNumNodesTotal(), 27);

  for (int i=0; i < 8; ++i)
    EXPECT_EQ(ref_el->getNodeIndex(0, i), i);

  for (int i=0; i < 12; ++i)
    EXPECT_EQ(ref_el->getNodeIndex(1, i, 0), i + 8);

  for (int i=0; i < 6; ++i)
    EXPECT_EQ(ref_el->getNodeIndex(2, i, 0), i + 20);

  EXPECT_EQ(ref_el->getNodeIndex(3, 0, 0), 26);
}


TEST(ReferenceElementHex, Cubic)
{
  auto ref_el = reference_element::getLagrangeHexReferenceElement(3);

  EXPECT_EQ(ref_el->getNumNodes(0), 1);
  EXPECT_EQ(ref_el->getNumNodes(1), 2);
  EXPECT_EQ(ref_el->getNumNodes(2), 4);
  EXPECT_EQ(ref_el->getNumNodes(3), 8);
  EXPECT_EQ(ref_el->getNumNodesTotal(), 64);

  int idx = 0;
  for (int i=0; i < 8; ++i)
    EXPECT_EQ(ref_el->getNodeIndex(0, i), idx++);

  for (int i=0; i < 12; ++i)
    for (int node=0; node < 2; ++node)
      EXPECT_EQ(ref_el->getNodeIndex(1, i, node), idx++);

  for (int i=0; i < 6; ++i)
    for (int node=0; node < 4; ++node)
      EXPECT_EQ(ref_el->getNodeIndex(2, i, node), idx++);


  for (int node=0; node < 8; ++node)
    EXPECT_EQ(ref_el->getNodeIndex(3, 0, node), idx++);
}
