
#include "ProjectDefs.h"
#include "lagrange_basis.h"
#include "mesh/mesh.h"
#include "mesh/reference_element.h"

namespace Heat {

// holds the basis values for the lagrange polynomials,
// presenting an interface as through there is no tensor-product
// structure
class BasisVals
{
  public: 
    BasisVals(const Mesh::TensorProductMapper& tp_mapper_in,
              const Mesh::TensorProductMapper& tp_mapper_out) :
    m_vals(lagrange_memoizer.getValues(tp_mapper_in.getXi(), tp_mapper_out.getXi())),
    m_derivs(lagrange_memoizer.getDerivs(tp_mapper_in.getXi(), tp_mapper_out.getXi())),
    m_rev_nodemap_in(boost::extents[tp_mapper_in.getNodemap().num_elements()][3]),
    m_rev_nodemap_out(boost::extents[tp_mapper_out.getNodemap().num_elements()][3])
    {
      getReverseNodemap(tp_mapper_in, m_rev_nodemap_in);
      getReverseNodemap(tp_mapper_out, m_rev_nodemap_out);
    }

    Real getValue(const int idx_in, const int idx_out)
    {
      int i_in = m_rev_nodemap_in[idx_in][0];  int i_out = m_rev_nodemap_out[idx_out][0];
      int j_in = m_rev_nodemap_in[idx_in][1];  int j_out = m_rev_nodemap_out[idx_out][1];
      int k_in = m_rev_nodemap_in[idx_in][2];  int k_out = m_rev_nodemap_out[idx_out][2];

      return m_vals[i_out][i_in] * m_vals[j_out][j_in] * m_vals[k_out][k_in];
    }

    void getDerivs(const int idx_in, const int idx_out, Real derivs[3])
    {
      int i_in = m_rev_nodemap_in[idx_in][0];  int i_out = m_rev_nodemap_out[idx_out][0];
      int j_in = m_rev_nodemap_in[idx_in][1];  int j_out = m_rev_nodemap_out[idx_out][1];
      int k_in = m_rev_nodemap_in[idx_in][2];  int k_out = m_rev_nodemap_out[idx_out][2];

      derivs[0] = m_derivs[i_out][i_in] * m_vals[j_out][j_in]   * m_vals[k_out][k_in];
      derivs[1] = m_vals[i_out][i_in]   * m_derivs[j_out][j_in] * m_vals[k_out][k_in];
      derivs[2] = m_vals[i_out][i_in]   * m_vals[j_out][j_in]   * m_derivs[k_out][k_in];
    }

    ArrayType<LocalIndex, 2>& getRevNodemapOut() { return m_rev_nodemap_out; }

  private:

    void getReverseNodemap(const Mesh::TensorProductMapper& mapper, ArrayType<LocalIndex, 2>& rev_nodemap)
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

    LagrangeMemoizer::RetType m_vals;
    LagrangeMemoizer::RetType m_derivs;
    ArrayType<LocalIndex, 2> m_rev_nodemap_in;  // maps linear index to tensor product triplet;
    ArrayType<LocalIndex, 2> m_rev_nodemap_out;
};


class BasisVals2D
{
  public:
    // face_nodemap is a nfaces per element x nnodes per face array giving the index of each volume node
    // that lies on the face
    BasisVals2D(const Mesh::TensorProductMapper& tp_mapper_in,
                std::vector<Real> pts_out,
                const ArrayType<LocalIndex, 2>& face_nodemap,
                Mesh::ReferenceElement* ref_el) :
      m_vals(boost::extents[face_nodemap.shape()[0]][face_nodemap.shape()[0]][pts_out.size()][pts_out.size()])
    {
      //TODO: there is a better way of doing this that requires less storage: both the input and the output are
      //      tensor product, so there should be a way to make m_vals 3 dimensional rather than 4 dimensional
      LagrangeBasis basis(tp_mapper_in.getXi());
      for (int face=0; face < getNumFaces(); ++face)
        for (int j=0; j < getNumNodesPerFace(); ++j)
          for (int ki=0; ki < getNumOutputNodes(); ++ki)
            for (int kj=0; kj < getNumOutputNodes(); ++kj)
            {
              std::array<Real, 2> xi_face{pts_out[ki], pts_out[kj]};
              std::array<Real, 3> xi_vol;
              ref_el->computeElementXi(face, xi_face.data(), xi_vol.data());


              auto tp_indices = getTPIndices(tp_mapper_in, face_nodemap[face][j]);
              m_vals[face][j][ki][kj] = basis.evalPoly(tp_indices[0], xi_vol[0]) * 
                                        basis.evalPoly(tp_indices[1], xi_vol[1]) *
                                        basis.evalPoly(tp_indices[2], xi_vol[2]);

            }
    }

    // face is the face of the element
    // face_node is the index of the solution node in the range [0, numSolNodesPerFace-1]
    // face_i and face_j are the tensor product indices of the face point
    Real getValue(int face, int face_node, int face_i, int face_j)
    {
      return m_vals[face][face_node][face_i][face_j];
    }

    int getNumFaces() const { return m_vals.shape()[0]; }

    int getNumNodesPerFace() const { return m_vals.shape()[1]; }

    int getNumOutputNodes() const { return m_vals.shape()[0]; }


  private:

    // get tensor product indicies of given volume node
    std::array<int, 3> getTPIndices(const Mesh::TensorProductMapper& tp_mapper_in, int node)
    {
      auto& nodemap = tp_mapper_in.getNodemap();
      for (int i=0; i < tp_mapper_in.getNumTPPoints(); ++i)
        for (int j=0; j < tp_mapper_in.getNumTPPoints(); ++j)
          for (int k=0; k < tp_mapper_in.getNumTPPoints(); ++k)
            if (nodemap[i][j][k] == node)
              return {i, j, k};

      throw std::runtime_error("unable to find specified volume node");
    }
    
    ArrayType<Real, 4> m_vals;  // nfaces per element x numNodesPerFace x num tensor product output nodes x num tensorProduct output nodes
    

};

}