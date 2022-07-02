#include "physics/heat/source_term.h"
#include "physics/heat/basis_vals.h"

namespace Heat {

void computeSourceTerm(const HeatEquation& physics, Real t, DiscVectorPtr rhs)
{
  //TODO: skip source term if it is zero/absent
  auto disc = physics.getDiscretization();

  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto vol_disc = disc->getVolDisc(i);
    auto src_term = physics.getSourceTerm(i);
    auto& rhs_arr = rhs->getArray(i);

    computeSourceTerm(vol_disc, src_term, t, rhs_arr);
  }

  rhs->markArrayModified();
}

void computeSourceTerm(const VolDiscPtr vol_disc, SourceTermPtr src, Real t,
                       ArrayType<Real, 2>& rhs_arr)
{
  auto& detJ  = vol_disc->detJ;
  BasisVals basis_vals(vol_disc->vol_group.getTPMapperSol(), Mesh::TensorProductMapper(vol_disc->quad.getPoints()));
  std::vector<Real> src_vals(vol_disc->getNumQuadPtsPerElement());

  auto& rev_nodemap = basis_vals.getRevNodemapOut();
  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    src->getValues(el, t, src_vals.data());

    for (int q=0; q < vol_disc->getNumQuadPtsPerElement(); ++q)
    {
      int k_i = rev_nodemap[q][0]; int k_j = rev_nodemap[q][1]; int k_k = rev_nodemap[q][2];
      Real weight = vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);

      for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      {
        rhs_arr[el][i] += basis_vals.getValue(i, q) * weight * src_vals[q] / detJ[el][q];
      }
    }
  }
}

}
