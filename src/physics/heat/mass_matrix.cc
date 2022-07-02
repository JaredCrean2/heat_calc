#include "physics/heat/mass_matrix.h"
#include "physics/heat/basis_vals.h"

namespace Heat {

void applyMassMatrix(const HeatEquation& physics, DiscVectorPtr vec_in, DiscVectorPtr vec_out)
{
  if (!vec_in->isArrayCurrent())
    vec_in->syncVectorToArray();

  auto disc = physics.getDiscretization();

  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto vol_disc = disc->getVolDisc(i);
    auto& params  = physics.getVolumeGroupParams(i);
    auto dofs     = disc->getDofNumbering();
    auto& arr_in  = vec_in->getArray(i);
    auto& arr_out = vec_out->getArray(i);

    applyMassMatrix(vol_disc, params, dofs, arr_in, arr_out);
  }

  vec_out->markArrayModified();
}

void applyMassMatrix(const VolDiscPtr vol_disc, const VolumeGroupParams& params,  const DofNumberingPtr dof_numbering,
                     const ArrayType<Real, 2>& arr_in, ArrayType<Real, 2>& arr_out)
{
  // for elements with dirichlet nodes, the form of the mass matrix is R^T M R,
  // where M is the full mass matrix, R is a diagonal matrix with 1 on the diagonal
  // for active dofs and 0 for dirichlet dofs.
  // Replacing arr_in with zeros for the non-active dofs is equivalent to R * arr_in,
  // and the transformation from array to vector form of arr_out is equivalent to R^T

  auto& detJ  = vol_disc->detJ;
  auto& tp_mapper_sol = vol_disc->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_disc->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);

  const auto& dof_nums = dof_numbering->getDofs(vol_disc);
  ArrayType<Real, 1> u_sol(boost::extents[vol_disc->getNumSolPtsPerElement()]);
  ArrayType<Real, 1> u_quad(boost::extents[vol_disc->getNumQuadPtsPerElement()]);
  auto& rev_nodemap = basis_vals.getRevNodemapOut();
  auto alpha = params.rho * params.Cp;

  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      u_sol[i] = dof_numbering->isDofActive(dof_nums[el][i]) ? arr_in[el][i] : 0;

    auto u_i = arr_in[boost::indices[el][range()]];
    vol_disc->interp_sq_flat_to_flat.interpolateVals(u_i, u_quad);

    for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
    {
      int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
      Real weight = vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);
      u_quad[k] = alpha * u_quad[k] * weight / detJ[el][k];
    }

    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
    {
      arr_out[el][i] = 0;
      for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
        arr_out[el][i] += basis_vals.getValue(i, k) * u_quad[k];
    }

  }
}

//-----------------------------------------------------------------------------
// computeMassMatrix
void computeMassMatrix(const HeatEquation& physics, linear_system::AssemblerPtr assembler)
{
  auto disc = physics.getDiscretization();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto vol_disc = disc->getVolDisc(i);
    auto& params  = physics.getVolumeGroupParams(i);

    computeMassMatrix(vol_disc, params, assembler);
  }

}

void computeMassMatrix(const VolDiscPtr vol_disc, const VolumeGroupParams& params, linear_system::AssemblerPtr assembler)
{
  auto& detJ  = vol_disc->detJ;
  Real rho_cp = params.rho * params.Cp;
  auto& tp_mapper_sol = vol_disc->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_disc->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);
  const auto& rev_nodemap = basis_vals.getRevNodemapOut();

  ArrayType<Real, 3> dN_dx(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumQuadPtsPerElement()][3]);
  ArrayType<Real, 2> dR_du(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumSolPtsPerElement()]);

  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      for (int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
      {
        dR_du[i][j] = 0;

        for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
        {
          Real Ni = basis_vals.getValue(i, k);
          Real Nj = basis_vals.getValue(j, k);

          int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
          //TODO: store 1/detJ in an array to avoid the division here
          Real weight = rho_cp * vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k) / detJ[el][k];
          dR_du[i][j] += Ni * weight * Nj;
        } 
      }

    assembler->assembleVolume(vol_disc->getIdx(), el, dR_du);
  }
}

}
