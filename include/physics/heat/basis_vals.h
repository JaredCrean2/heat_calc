
#include "ProjectDefs.h"
#include "lagrange_basis.h"
#include "mesh/mesh.h"
#include "mesh/reference_element.h"
#include "mesh/reference_element_geometry.h"

namespace Heat {

// holds the basis values for the lagrange polynomials,
// presenting an interface as through there is no tensor-product
// structure
class BasisVals
{
  public: 
    BasisVals(const Mesh::TensorProductMapper& tp_mapper_in,
              const Mesh::TensorProductMapper& tp_mapper_out);

    Real getValue(const int idx_in, const int idx_out) const
    {
      int i_in = m_rev_nodemap_in[idx_in][0];  int i_out = m_rev_nodemap_out[idx_out][0];
      int j_in = m_rev_nodemap_in[idx_in][1];  int j_out = m_rev_nodemap_out[idx_out][1];
      int k_in = m_rev_nodemap_in[idx_in][2];  int k_out = m_rev_nodemap_out[idx_out][2];

      return m_vals[i_out][i_in] * m_vals[j_out][j_in] * m_vals[k_out][k_in];
    }

    void getDerivs(const int idx_in, const int idx_out, Real derivs[3]) const
    {
      int i_in = m_rev_nodemap_in[idx_in][0];  int i_out = m_rev_nodemap_out[idx_out][0];
      int j_in = m_rev_nodemap_in[idx_in][1];  int j_out = m_rev_nodemap_out[idx_out][1];
      int k_in = m_rev_nodemap_in[idx_in][2];  int k_out = m_rev_nodemap_out[idx_out][2];

      derivs[0] = m_derivs[i_out][i_in] * m_vals[j_out][j_in]   * m_vals[k_out][k_in];
      derivs[1] = m_vals[i_out][i_in]   * m_derivs[j_out][j_in] * m_vals[k_out][k_in];
      derivs[2] = m_vals[i_out][i_in]   * m_vals[j_out][j_in]   * m_derivs[k_out][k_in];
    }

    const ArrayType<LocalIndex, 2>& getRevNodemapIn() const { return m_rev_nodemap_in; }

    const ArrayType<LocalIndex, 2>& getRevNodemapOut() const { return m_rev_nodemap_out; }

  private:

    void getReverseNodemap(const Mesh::TensorProductMapper& mapper, ArrayType<LocalIndex, 2>& rev_nodemap);

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
                Mesh::REPtr ref_el);

    // face is the face of the element
    // face_node is the index of the solution node in the range [0, numSolNodesPerFace-1]
    // face_i and face_j are the tensor product indices of the face point
    Real getValue(int face, int face_node, int face_i, int face_j)
    {
      return m_vals[face][face_node][face_i][face_j];
    }

    int getNumFaces() const { return m_vals.shape()[0]; }

    int getNumNodesPerFace() const { return m_vals.shape()[1]; }

    int getNumOutputNodes() const { return m_vals.shape()[2]; }


  private:

    // get tensor product indicies of given volume node
    std::array<int, 3> getTPIndices(const Mesh::TensorProductMapper& tp_mapper_in, int node);
    
    ArrayType<Real, 4> m_vals;  // nfaces per element x numNodesPerFace x num tensor product output nodes x num tensorProduct output nodes
};

}