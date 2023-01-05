#include "physics/heat/neumann_bc.h"
#include "physics/heat/basis_vals.h"
#include "physics/heat/helper_funcs.h"

namespace Heat {

void computeNeumannBC(const HeatEquation& physics, const Real t, DiscVectorPtr u, DiscVectorPtr rhs)
{
  const auto& neumann_bcs = physics.getNeumannBCs();
  for (auto& bc : neumann_bcs)
  {
    //std::cout << "doing BC " << i << std::endl;
    computeNeumannBC(bc, u, t, rhs);
  }

  rhs->markArrayModified();
}

//TODO: a face could be a ghost face.  Need weighting factor?
//      or maybe it is ok beause we throw away any non-owned nodes,
//      and the the owned nodes have the correct parallel value

void computeNeumannBC(NeumannBCPtr bc, DiscVectorPtr u, const Real t, DiscVectorPtr rhs)
{
  // Note: we currently don't require that a given surface has a single volume as
  //       as its upward adjacency, so we have to look up a different volume disc
  //       for every face
  
  auto surf = bc->getSurfDisc();
  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::vector<Real> flux_vals(surf->getNumQuadPtsPerFace() * 3);
  Quadrature& quad = surf->quad;
  BasisVals2D basis(surf->face_group.getTPMapperSol(), quad.getPoints(), surf->face_group.getFaceNodesSol(), surf->face_group.ref_el_sol);

  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    //std::cout << "face = " << face << std::endl;
    auto& face_spec = surf->face_group.faces[face];
    auto& u_arr = u->getArray(face_spec.vol_group);
    auto& res_arr = rhs->getArray(face_spec.vol_group);

    auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];
    auto res_el = res_arr[boost::indices[face_spec.el_group][range()]];

    surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);

    bc->getValue(face, t, u_quad.data(), flux_vals.data());
    //for (int k=0; k < surf->getNumQuadPtsPerFace(); ++k)
    //  for (int d=0; d < 3; ++d)
    //    std::cout << "k = " << k << ", d = " << d << ", flux = " << flux_vals[k + surf->getNumQuadPtsPerFace()*d] << std::endl;

    for (int ki=0; ki < quad.getNumPoints(); ++ki)
      for (int kj=0; kj < quad.getNumPoints(); ++kj)
      {
        int k    = surf->quad_tp_nodemap[ki][kj];
        Real weight = quad.getWeight(ki) * quad.getWeight(kj);
        Real flux_normal = 0;
        for (int d=0; d < 3; ++d)
          flux_normal += surf->normals[face][k][d] * flux_vals[k + surf->getNumQuadPtsPerFace() * d];

        //std::cout << "k = " << k << ", flux_normal = " << flux_normal << std::endl;


        Real val = weight * flux_normal;
        for (int i=0; i < surf->getNumSolPtsPerFace(); ++i)
        {
          int node_sol = surf->face_group.getFaceNodesSol()[face_spec.face][i];
          res_arr[face_spec.el_group][node_sol] += basis.getValue(face_spec.face, i, ki, kj) * val;
        }          
      }
  }
}

void computeNeumannBCJacobian(const HeatEquation& physics, DiscVectorPtr u, Real t, linear_system::AssemblerPtr assembler)
{
  const auto& neumann_bcs = physics.getNeumannBCs();
  for (auto& bc : neumann_bcs)
  {
    if (bc->isNonlinear())
      computeNeumannBCJacobian(bc, u, t, assembler);
  }
}

void computeNeumannBCJacobian(NeumannBCPtr bc, DiscVectorPtr u, Real t, linear_system::AssemblerPtr assembler)
{
  // Note: we currently don't require that a given surface has a single volume as
  //       as its upward adjacency, so we have to look up a different volume disc
  //       for every face
  
  auto surf = bc->getSurfDisc();
  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::vector<Real> flux_vals_deriv(surf->getNumQuadPtsPerFace() * 3);
  ArrayType<Real, 2> dRdu(boost::extents[surf->getNumSolPtsPerFace()][surf->getNumSolPtsPerFace()]);
  Quadrature& quad = surf->quad;
  BasisVals2D basis(surf->face_group.getTPMapperSol(), quad.getPoints(), surf->face_group.getFaceNodesSol(), surf->face_group.ref_el_sol);

  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    auto& face_spec = surf->face_group.faces[face];
    auto& u_arr = u->getArray(face_spec.vol_group);

    auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];
    //auto res_el = res_arr[boost::indices[face_spec.el_group][range()]];
    zeroMatrix(dRdu);

    surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);
    bc->getValueDeriv(face, t, u_quad.data(), flux_vals_deriv.data());

    //TODO: rewrite this for vectorization
    for (int ki=0; ki < quad.getNumPoints(); ++ki)
      for (int kj=0; kj < quad.getNumPoints(); ++kj)
      {
        int k       = surf->quad_tp_nodemap[ki][kj];
        Real weight = quad.getWeight(ki) * quad.getWeight(kj);
        Real flux_normal = 0;
        for (int d=0; d < 3; ++d)
          flux_normal += surf->normals[face][k][d] * flux_vals_deriv[k + surf->getNumQuadPtsPerFace() * d];

        Real val = weight * flux_normal;

        // TODO: it would be better if kj was the innermost loop
        for (int i=0; i < surf->getNumSolPtsPerFace(); ++i)
          for (int j=0; j < surf->getNumSolPtsPerFace(); ++j)
            dRdu[i][j] += basis.getValue(face_spec.face, i, ki, kj) * val * basis.getValue(face_spec.face, j, ki, kj);
      }

    assembler->assembleFace(surf->getIdx(), face, dRdu);
  }
}

}
