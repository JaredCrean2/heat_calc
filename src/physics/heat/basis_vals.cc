#include "physics/heat/basis_vals.h"

namespace Heat {

BasisVals::BasisVals(const Mesh::TensorProductMapper& tp_mapper_in,
                     const Mesh::TensorProductMapper& tp_mapper_out) :
  m_vals(lagrange_memoizer.getValues(tp_mapper_in.getXi(), tp_mapper_out.getXi())),
  m_derivs(lagrange_memoizer.getDerivs(tp_mapper_in.getXi(), tp_mapper_out.getXi())),
  m_rev_nodemap_in(boost::extents[tp_mapper_in.getNodemap().num_elements()][3]),
  m_rev_nodemap_out(boost::extents[tp_mapper_out.getNodemap().num_elements()][3])
{
  getReverseNodemap(tp_mapper_in, m_rev_nodemap_in);
  getReverseNodemap(tp_mapper_out, m_rev_nodemap_out);
}


void BasisVals::getReverseNodemap(const Mesh::TensorProductMapper& mapper, ArrayType<LocalIndex, 2>& rev_nodemap)
{
  auto& nodemap = mapper.getNodemap();
  for (unsigned int i=0; i < nodemap.shape()[0]; ++i)
    for (unsigned int j=0; j < nodemap.shape()[1]; ++j)
      for (unsigned int k=0; k < nodemap.shape()[2]; ++k)
      {
        rev_nodemap[nodemap[i][j][k]][0] = i;
        rev_nodemap[nodemap[i][j][k]][1] = j;
        rev_nodemap[nodemap[i][j][k]][2] = k;
      }
}


BasisVals2D::BasisVals2D(const Mesh::TensorProductMapper& tp_mapper_in,
            std::vector<Real> pts_out,
            const ArrayType<LocalIndex, 2>& face_nodemap,
            Mesh::REPtr ref_el) :
  m_vals(boost::extents[face_nodemap.shape()[0]][face_nodemap.shape()[1]][pts_out.size()][pts_out.size()])
{
  //TODO: there is a better way of doing this that requires less storage: both the input and the output are
  //      tensor product, so there should be a way to make m_vals 3 dimensional rather than 4 dimensional
  LagrangeBasis basis(tp_mapper_in.getXi());
  for (int face=0; face < getNumFaces(); ++face)
  {
    auto ge_el = ref_el->getElement();
    auto ge_face = ref_el->getFace(face);
    for (int j=0; j < getNumNodesPerFace(); ++j)
      for (int ki=0; ki < getNumOutputNodes(); ++ki)
        for (int kj=0; kj < getNumOutputNodes(); ++kj)
        {
          std::array<Real, 3> xi_face{pts_out[ki], pts_out[kj], 0};
          //TODO: this is likely a performance problem
          auto xi_vol = reference_element::reclassifyPoint(ge_face, xi_face, ge_el);
          //ref_el->computeElementXi(face, xi_face.data(), xi_vol.data());

          auto tp_indices = getTPIndices(tp_mapper_in, face_nodemap[face][j]);
          m_vals[face][j][ki][kj] = basis.evalPoly(tp_indices[0], xi_vol[0]) * 
                                    basis.evalPoly(tp_indices[1], xi_vol[1]) *
                                    basis.evalPoly(tp_indices[2], xi_vol[2]);
    }
  }
}

// get tensor product indicies of given volume node
std::array<int, 3> BasisVals2D::getTPIndices(const Mesh::TensorProductMapper& tp_mapper_in, int node)
{
  auto& nodemap = tp_mapper_in.getNodemap();
  for (int i=0; i < tp_mapper_in.getNumTPPoints(); ++i)
    for (int j=0; j < tp_mapper_in.getNumTPPoints(); ++j)
      for (int k=0; k < tp_mapper_in.getNumTPPoints(); ++k)
        if (nodemap[i][j][k] == node)
          return {i, j, k};

  throw std::runtime_error("unable to find specified volume node");
}


} // namespace