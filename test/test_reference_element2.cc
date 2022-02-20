#include "gtest/gtest.h"
#include "mesh/reference_element_interface.h"
#include "mesh/reference_element_geometry_interface.h"

namespace {

using REPtr = std::shared_ptr<reference_element::ReferenceElement>;

void testTensorProductXi(REPtr ref_el, const std::vector<Real>& tp_xi_exact)
{
  auto& tp_xi = ref_el->getTensorProductXi();

  EXPECT_EQ(tp_xi.size(), tp_xi_exact.size());
  for (size_t i=0; i < tp_xi_exact.size(); ++i)
    EXPECT_NEAR(tp_xi[i], tp_xi_exact[i], 1e-13);
}

void testTensorProductNodemap(REPtr ref_el)
{
  std::vector<int> node_indices;
  auto& tp_nodemap = ref_el->getTPNodemap();

  for (int i=0; i < ref_el->getNumNodesTP(); ++i)
    for (int j=0; j < ref_el->getNumNodesTP(); ++j)
      for (int k=0; k < ref_el->getNumNodesTP(); ++k)
        node_indices.push_back(tp_nodemap[i][j][k]);

  std::sort(node_indices.begin(), node_indices.end());
  EXPECT_EQ(node_indices.size(), ref_el->getNumNodesTotal());
  EXPECT_EQ(node_indices.front(), 0);
  for (unsigned int i=1; i < node_indices.size(); ++i)
    EXPECT_EQ(node_indices[i], node_indices[i-1] + 1);
}

void testNodeXi(REPtr ref_el, const std::vector<Real>& tp_xi_exact)
{
  auto& tp_nodemap = ref_el->getTPNodemap();
  auto& node_xi    = ref_el->getNodeXi();
  
  for (int i=0; i < ref_el->getNumNodesTP(); ++i)
    for (int j=0; j < ref_el->getNumNodesTP(); ++j)
      for (int k=0; k < ref_el->getNumNodesTP(); ++k)
      {
        int node_idx = tp_nodemap[i][j][k];
        EXPECT_NEAR(node_xi[node_idx][0], tp_xi_exact[i], 1e-13);
        EXPECT_NEAR(node_xi[node_idx][1], tp_xi_exact[j], 1e-13);
        EXPECT_NEAR(node_xi[node_idx][2], tp_xi_exact[k], 1e-13);
      }
}

void testXiRange(REPtr ref_el)
{
  auto p = ref_el->getXiRange();

  EXPECT_NEAR(p.first,  0, 1e-13);
  EXPECT_NEAR(p.second, 1, 1e-13);
}

template <typename T1, typename T2>
bool compareArrays(const T1& a, const T2& b, Real tol)
{
  return std::abs(a[0] - b[0]) < tol &&
         std::abs(a[1] - b[1]) < tol &&
         std::abs(a[2] - b[2]) < tol;
}

void testNormals(REPtr ref_el)
{
  std::vector< std::array<Real, 3> > normals_ex{ {0, 0, -1}, {0, -1, 0}, {1, 0, 0},
                                                 {0, 1, 0},  {-1, 0, 0}, {0, 0, 1}};
  auto& normals = ref_el->getNormals();
  EXPECT_EQ(normals_ex.size(), ref_el->getNumFaces());
  for (int i=0; i < ref_el->getNumFaces(); ++i)
  {
    auto normals_i = normals[boost::indices[i][range()]];
    EXPECT_TRUE(compareArrays(normals_i, normals_ex[i], 1e-13));
  }
}

}


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

  std::vector<Real> tp_xi_exact{0.0, 1.0};
  testTensorProductXi(ref_el, tp_xi_exact);
  testTensorProductNodemap(ref_el);
  testNodeXi(ref_el, tp_xi_exact);
  testXiRange(ref_el);
  testNormals(ref_el);

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

  std::vector<Real> tp_xi_exact{0.0, 0.5, 1.0};
  testTensorProductXi(ref_el, tp_xi_exact);
  testTensorProductNodemap(ref_el);
  testNodeXi(ref_el, tp_xi_exact);
  testXiRange(ref_el);
  testNormals(ref_el);

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

  std::vector<Real> tp_xi_exact{0.0, 1.0/3.0, 2.0/3.0, 1.0};
  testTensorProductXi(ref_el, tp_xi_exact);
  testTensorProductNodemap(ref_el);
  testNodeXi(ref_el, tp_xi_exact);
  testXiRange(ref_el);
  testNormals(ref_el);

}
