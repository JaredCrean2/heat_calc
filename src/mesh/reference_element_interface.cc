#include "mesh/reference_element_interface.h"
#include "mesh/reference_element_nodes.h"
#include "mesh/reference_element_geometry_hex.h"

namespace reference_element {

ReferenceElement::ReferenceElement(const ReferenceElementGeometry& geom, const ReferenceNodes& nodes) :
  m_geom(geom),
  m_nodes(nodes)
{
  auto tp_nodemap = computeTPNodemap(*this);
  auto normals    = computeNormals(*this);

  m_tp_nodemap.resize(boost::extents[tp_nodemap.shape()[0]][tp_nodemap.shape()[1]][tp_nodemap.shape()[2]]);
  m_tp_nodemap = tp_nodemap;

  m_normals.resize(boost::extents[normals.shape()[0]][normals.shape()[1]]);
  m_normals = normals;
}

std::vector<Real> get1DXi(ReferenceElement& ref_el)
{
  int ref_edge = 0;
  auto verts = ref_el.getDef().edge_verts[ref_edge];
  auto vert1 = ref_el.getVert(verts[0]);
  auto vert2 = ref_el.getVert(verts[1]);
  auto edge  = ref_el.getEdge(ref_edge);

  std::vector<Real> xi_vals;
  auto pt = edge->computeFromDownCoordinates(vert1, Point{0, 0, 0});
  xi_vals.push_back(pt[0]);

  for (int i=0; i < ref_el.getNumNodes(1); ++i)
    xi_vals.push_back(ref_el.getNode(1, i).xi[0]);

  pt = edge->computeFromDownCoordinates(vert2, Point{0, 0, 0});
  xi_vals.push_back(pt[0]);
  assertAlways(std::is_sorted(xi_vals.begin(), xi_vals.end()), "edge is not oriented from vert1 to vert2");

  return xi_vals;
}

int searchNode(ReferenceElement& ref_el, const Point& pt)
{
  GEPtr el = ref_el.getElement();
  for (int dim=0; dim <= 3; ++dim)
    for (int i=0; i < ref_el.getNumEntities(dim); ++i)
    {
      GEPtr entity = ref_el.getEntity(dim, i);
      for (int node=0; node < ref_el.getNumNodes(dim); ++node)
      {
        auto node_obj = ref_el.getNode(dim, node);
        auto pt_node = reclassifyPoint(entity, node_obj.xi, el);

        if (std::abs(pt[0] - pt_node[0]) < 1e-13 && std::abs(pt[1] - pt_node[1]) < 1e-13 &&
            std::abs(pt[2] - pt_node[2]) < 1e-13)
          return ref_el.getNodeIndex(dim, i, node);

      }
    }

  throw std::invalid_argument("could not find point");
}


ArrayType<LocalIndex, 3> computeTPNodemap(ReferenceElement& ref_el)
{
  auto xi_vals = get1DXi(ref_el);
  int nnodes = xi_vals.size();

  ArrayType<LocalIndex, 3> nodemap(boost::extents[nnodes][nnodes][nnodes]);
  for (int k=0; k < nnodes; ++k)
    for (int j=0; j < nnodes; ++j)
      for (int i=0; i < nnodes; ++i)
      {
        Point pt{xi_vals[i], xi_vals[j], xi_vals[k]};
        nodemap[i][j][k] = searchNode(ref_el, pt);
      }

  return nodemap;
}


ArrayType<Real, 2> computeNormals(ReferenceElement& ref_el)
{
  int nfaces = ref_el.getNumEntities(2);
  ArrayType<Real, 2> normals(boost::extents[nfaces][3]);
  auto el = ref_el.getElement();

  for (int i=0; i < nfaces; ++i)
  {
    auto face = ref_el.getFace(i);
    auto verts = getDownward(face, 0);

    auto node = ref_el.getNode(0, 0);

    auto pt0 = reclassifyPoint(verts[0], node.xi, el);
    auto pt1 = reclassifyPoint(verts[1], node.xi, el);
    auto pt3 = reclassifyPoint(verts[3], node.xi, el);

    auto b1 = pt1 - pt0;
    auto b2 = pt3 - pt0;
    auto norm = cross(b1, b2);

    // make sure normal is oriented outwards by dotting normal with
    // vector from pt0 to a vertex not on the face

    GEPtr other_vert;
    for (int j=0; j < ref_el.getNumVerts(); ++j)
      if (std::find(verts.begin(), verts.end(), ref_el.getVert(j)) == verts.end())
      {
        other_vert = ref_el.getVert(j);
        break;
      }

    auto pt_other = reclassifyPoint(other_vert, node.xi, el);
    auto b3 = pt_other - pt0;

    auto val = dot(norm, b3);
    assert(std::abs(val) > 1e-7);
    if (val < 0)
      norm = norm * -1;

    for (int j=0; j < 3; ++j)
      normals[i][j] = norm[j];
  }

  return normals;
}


std::shared_ptr<ReferenceElement> getLagrangeHexReferenceElement(int degree)
{
  auto def = getStandardReferenceElementDef();
  ReferenceElementGeometry ref_geom(def);
  auto nodes = getLagrangeHexNodes(degree);

  return std::make_shared<ReferenceElement>(ref_geom, nodes);
}

}