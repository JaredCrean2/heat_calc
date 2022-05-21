#include "gtest/gtest.h"
#include "test_helper.h"
#include "mesh/reference_element_geometry.h"
#include "mesh/reference_element_geometry_hex.h"
#include "mesh/reference_element_geometry_interface.h"

using Point = reference_element::Point;

bool isNear(const Point& pt1, const Point& pt2, Real tol)
{
  return std::abs(pt1[0] - pt2[0]) < tol &&
         std::abs(pt1[1] - pt2[1]) < tol &&
         std::abs(pt1[0] - pt2[0]) < tol;
}

TEST(ReferenceGeometry, Verifier)
{
  SERIAL_ONLY();
  EXPECT_NO_THROW(reference_element::checkReferenceElementDef(reference_element::getStandardReferenceElementDef()));
}


TEST(ReferenceGeometry, Counts)
{
  SERIAL_ONLY();

  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  EXPECT_EQ(ref_el.getNumVerts(), 8);
  EXPECT_EQ(ref_el.getNumEdges(), 12);
  EXPECT_EQ(ref_el.getNumFaces(), 6);

  for (int i=0; i < ref_el.getNumVerts(); ++i)
  {
    EXPECT_NE(ref_el.getVert(i), nullptr);
    EXPECT_EQ(ref_el.getVert(i)->getDimension(), 0);
    EXPECT_EQ(ref_el.getVert(i)->getId(), i);
  }

  for (int i=0; i < ref_el.getNumEdges(); ++i)
  {
    EXPECT_NE(ref_el.getEdge(i), nullptr);
    EXPECT_EQ(ref_el.getEdge(i)->getDimension(), 1);
    EXPECT_EQ(ref_el.getEdge(i)->getId(), i);
  }

  for (int i=0; i < ref_el.getNumFaces(); ++i)
  {
    EXPECT_NE(ref_el.getFace(i), nullptr);
    EXPECT_EQ(ref_el.getFace(i)->getDimension(), 2);
    EXPECT_EQ(ref_el.getFace(i)->getId(), i);
  }

  EXPECT_NE(ref_el.getElement(), nullptr);
  EXPECT_EQ(ref_el.getElement()->getDimension(), 3);
  EXPECT_EQ(ref_el.getElement()->getId(), 0);
}


TEST(ReferenceGeometry, VertCoords)
{
  SERIAL_ONLY();

  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  Point pt_in{0, 0, 0};
  std::vector<Point> pt_exact{ {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

  auto el = ref_el.getElement();
  for (int i=0; i < ref_el.getNumVerts(); ++i)
  {
    auto pt_out = reference_element::reclassifyPoint(ref_el.getVert(i), pt_in, el);
    EXPECT_TRUE(isNear(pt_out, pt_exact[i], 1e-13));
  }
}


TEST(ReferenceGeometry, EdgeCoords)
{
  SERIAL_ONLY();

  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  Point pt_in{0.25, 0, 0};
  std::vector<Point> pt_exact{ {0.25, 0, 0}, {1, 0.25, 0}, {0.75, 1, 0}, {0, 0.25, 0},
                               {0.25, 0, 1}, {1, 0.25, 1}, {0.75, 1, 1}, {0, 0.25, 1},
                               {0, 0, 0.25}, {1, 0, 0.25}, {1, 1, 0.25}, {0, 1, 0.25}};

  auto el = ref_el.getElement();
  for (int i=0; i < ref_el.getNumEdges(); ++i)
  {
    auto pt_out = reference_element::reclassifyPoint(ref_el.getEdge(i), pt_in, el);
    EXPECT_TRUE(isNear(pt_out, pt_exact[i], 1e-13));
  }
}


TEST(ReferenceGeometry, FaceCoords)
{
  SERIAL_ONLY();

  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  Point pt_in{0.25, 0.35, 0};
  std::vector<Point> pt_exact{ {0.25, 0.35, 0}, {0.25, 0, 0.35}, {1, 0.25, 0.35},
                               {0.75, 1, 0.35}, {0, 0.25, 0.35}, {0.25, 0.35, 1}};
;

  auto el = ref_el.getElement();
  for (int i=0; i < ref_el.getNumFaces(); ++i)
  {
    auto pt_out = reference_element::reclassifyPoint(ref_el.getFace(i), pt_in, el);
    EXPECT_TRUE(isNear(pt_out, pt_exact[i], 1e-13));
  }
}


TEST(ReferenceGeometry, getDownwardElToFace)
{
  SERIAL_ONLY();
  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  auto el = ref_el.getElement();

  auto faces = reference_element::getDownward(el, 2);
  EXPECT_EQ(faces.size(), static_cast<size_t>(6));
  for (int i=0; i < 6; ++i)
  {
    EXPECT_EQ(faces[i]->getDimension(), 2);
    EXPECT_EQ(faces[i]->getId(), i);
  }
}

TEST(ReferenceGeometry, getDownwardFaceToEdge)
{
  SERIAL_ONLY();

  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  for (int i=0; i < 6; ++i)
  {
    auto face = ref_el.getFace(i);

    auto edges = reference_element::getDownward(face, 1);
    auto& faces_def = ref_el.getDef().face_edges;
    EXPECT_EQ(edges.size(), static_cast<size_t>(4));
    for (int j=0; j < 4; ++j)
    {
      EXPECT_EQ(edges[j]->getDimension(), 1);
      EXPECT_EQ(edges[j]->getId(), faces_def[i][j]);
    }
  }
}


TEST(ReferenceGeometry, getDownwardElToEdge)
{
  SERIAL_ONLY();

  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  std::vector<int> edges_expected = {0, 1, 2, 3, 9, 4, 8, 10, 5, 11, 6, 7};

  auto el = ref_el.getElement();
  auto edges = reference_element::getDownward(el, 1);
  EXPECT_EQ(edges.size(), static_cast<size_t>(12));
  for (int i=0; i < 12; ++i)
  {
    EXPECT_EQ(edges[i]->getDimension(), 1);
    EXPECT_EQ(edges[i]->getId(), edges_expected[i]);
  }
}


TEST(ReferenceGeometry, getDownwardElToVert)
{
  SERIAL_ONLY();
  
  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  std::vector<int> verts_expected = {0, 1, 2, 3, 5, 4, 6, 7};

  auto el = ref_el.getElement();
  auto verts = reference_element::getDownward(el, 0);
  EXPECT_EQ(verts.size(), static_cast<size_t>(8));
  for (int i=0; i < 8; ++i)
  {
    EXPECT_EQ(verts[i]->getDimension(), 0);
    EXPECT_EQ(verts[i]->getId(), verts_expected[i]);
  }
}


