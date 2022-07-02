#include "physics/heat/dirichlet_bc.h"
#include "physics/heat/basis_vals.h"

namespace Heat {

void applyDirichletValues(const HeatEquation& physics, const Real t, DiscVectorPtr u)
{
  for (auto bc : physics.getDirichletBCs())
    applyDirichletValues(bc, t, u);
}

void computeUnsteadyDirichletBC(const HeatEquation& physics, const Real t, DiscVectorPtr rhs)
{
  auto disc = physics.getDiscretization();
  std::vector<VolDiscPtr> vol_discs(disc->getNumVolDiscs());
  std::vector<const VolumeGroupParams*> params(disc->getNumVolDiscs());
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    vol_discs[i] = disc->getVolDisc(i);
    params[i] = &(physics.getVolumeGroupParams(i));
  }

  for (auto& bc : physics.getDirichletBCs())
    if (bc->getIsUnsteady())
      computeUnsteadyDirichletBC(bc, vol_discs, params, t, rhs);
}


void computeUnsteadyDirichletBC(DirichletBCPtr bc, const std::vector<VolDiscPtr>& vol_discs,
                                const std::vector<const VolumeGroupParams*>& vol_group_params,
                                const Real t, DiscVectorPtr rhs)
{
  // Note: we currently don't require that a given surface has a single volume as
  //       as its upward adjacency, so we have to look up a different volume disc
  //       for every face
  
  auto surf = bc->getSurfDisc();
  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::vector<Real> dudt_vals_face(surf->getNumSolPtsPerFace());
  auto& face_to_volume_nodes = surf->face_group.ref_el_sol->getFaceNodes();
  auto& tp_mapper_sol = vol_discs[0]->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_discs[0]->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);
  const auto& rev_nodemap = basis_vals.getRevNodemapOut();

  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    auto& face_spec = surf->face_group.faces[face];
    auto vol_disc = vol_discs[face_spec.vol_group];
    int el = face_spec.el_group;
    Real rho_cp = vol_group_params[face_spec.vol_group]->rho * vol_group_params[face_spec.vol_group]->Cp;
    auto& detJ = vol_disc->detJ;
    auto& res_arr = rhs->getArray(face_spec.vol_group);

    bc->getValueDt(face, t, dudt_vals_face.data());

    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)  //TODO: could skip Dirichlet nodes
      for (int face_pt=0; face_pt < face_to_volume_nodes.shape()[1]; ++face_pt)
      {
        int j = face_to_volume_nodes[face_spec.face][face_pt];
        Real dudt_val = dudt_vals_face[face_pt];
        for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
        {
          Real Ni = basis_vals.getValue(i, k);
          Real Nj = basis_vals.getValue(j, k);

          int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
          //TODO: store 1/detJ in an array to avoid the division here
          Real weight = rho_cp * vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k) / detJ[el][k];
          res_arr[el][i] -= Ni * weight * Nj * dudt_val;     
        }
      }
  }
}




}
