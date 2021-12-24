
#include "ProjectDefs.h"
#include "lagrange_basis.h"
#include "mesh/mesh.h"

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

}