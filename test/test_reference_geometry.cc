#include "gtest/gtest.h"
#include "mesh/reference_element_geometry.h"
#include "mesh/reference_element_geometry_types.h"

using Point = reference_element::Point;

bool isNear(const Point& pt1, const Point& pt2, Real tol)
{
  return std::abs(pt1[0] - pt2[0]) < tol &&
         std::abs(pt1[1] - pt2[1]) < tol &&
         std::abs(pt1[0] - pt2[0]) < tol;
}

TEST(ReferenceGeometry, Verifier)
{
  EXPECT_NO_THROW(reference_element::checkReferenceElementDef(reference_element::getStandardReferenceElementDef()));
}


TEST(ReferenceGeometry, Counts)
{
  reference_element::ReferenceElementGeometry ref_el(reference_element::getStandardReferenceElementDef());

  EXPECT_EQ(ref_el.getNumVerts(), 8);
  EXPECT_EQ(ref_el.getNumEdges(), 12);
  EXPECT_EQ(ref_el.getNumFaces(), 6);

  for (int i=0; i < ref_el.getNumVerts(); ++i)
    EXPECT_NE(ref_el.getVert(i), nullptr);

  for (int i=0; i < ref_el.getNumEdges(); ++i)
      EXPECT_NE(ref_el.getEdge(i), nullptr);

  for (int i=0; i < ref_el.getNumFaces(); ++i)
      EXPECT_NE(ref_el.getFace(i), nullptr);

  EXPECT_NE(ref_el.getElement(), nullptr);
}


TEST(ReferenceGeometry, VertCoords)
{
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