#include "physics/heat/HeatEquation.h"
#include "discretization/disc_vector.h"
#include "physics/heat/basis_vals.h"
#include "mesh/mesh.h"

namespace Heat {


void HeatEquation::computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs)
{
  rhs->set(0);
  if (!(u->isArrayCurrent()))
    u->syncVectorToArray();

  std::cout << "before setting dirichlet values, u = " << std::endl;
  printArray(u);

  // apply Dirichlet values to array
  applyDirichletValues(*this, t, u);
  std::cout << "after setting dirichlet values, u = " << std::endl;
  printArray(u);

  // compute volume terms
  computeVolumeTerm(*this, u, rhs);

  std::cout << "\nafter volume term " << std::endl;
  printVector(rhs);

  // compute Neumann BC terms
  computeNeumannBC(*this, t, u, rhs);

  //std::cout << "\nafter Neumann term " << std::endl;
  //printArray(rhs);

  // compute source term
  computeSourceTerm(*this, t, rhs);

  std::cout << "\nafter source term " << std::endl;
  //printArray(rhs);
  printVector(rhs);
}

void HeatEquation::checkInitialization()
{
  PhysicsModel::checkInitialization();
  
  if (m_params.size() != getDiscretization()->getNumVolDiscs())
    throw std::runtime_error("Incorrect number of Heat::VolumeGroupParams.  Should be equal to number of volume groups");
}


void applyDirichletValues(const HeatEquation& physics, const Real t, DiscVectorPtr u)
{
  for (auto bc : physics.getDirichletBCs())
    applyDirichletValues(bc, t, u);
}

//-----------------------------------------------------------------------------
// Volume term

void computeVolumeTerm(const HeatEquation& physics, DiscVectorPtr u, DiscVectorPtr rhs)
{
  auto disc = physics.getDiscretization();

  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto vol_disc = disc->getVolDisc(i);
    auto& params  = physics.getVolumeGroupParams(i);
    auto& u_arr   = u->getArray(i);
    auto& rhs_arr = rhs->getArray(i);

    computeVolumeTerm(vol_disc, params, u_arr, rhs_arr);
  }

  rhs->markArrayModified();
}


void computeVolumeTerm(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                       const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr)
{
  auto& dxidx = vol_disc->dxidx;
  auto& detJ  = vol_disc->detJ;
  Real alpha = params.kappa / (params.rho * params.Cp);
  BasisVals basis_vals(vol_disc->vol_group.getTPMapperSol(), Mesh::TensorProductMapper(vol_disc->quad.getPoints()));
  Real dNi_dxi[3], dNi_dx[3], dNj_dxi[3], dNj_dx[3];

  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
    {
      for (int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
        for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
        {
          // evaluate shape function derivatives dN/dxi at quadrature point
          basis_vals.getDerivs(i, k, dNi_dxi);
          basis_vals.getDerivs(j, k, dNj_dxi);


          // use metrics to rotate dN/dxi to dN/dx
          for (int d1=0; d1 < 3; ++d1)
          {
            dNi_dx[d1] = 0;
            dNj_dx[d1] = 0;
            for (int d2=0; d2 < 3; ++d2)
            {
              dNi_dx[d1] += dxidx[el][k][d2][d1] * dNi_dxi[d2];
              dNj_dx[d1] += dxidx[el][k][d2][d1] * dNj_dxi[d2];
            }
          }

          // compute dN_i/dx * weights * dN_j/dx * |J| * T_j
          auto& rev_nodemap = basis_vals.getRevNodemapOut();
          int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
          Real weight = vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);
          for (int d=0; d < 3; ++d)
          {
            rhs_arr[el][i] += alpha * dNi_dx[d] * weight * dNj_dx[d] * u_arr[el][j] / detJ[el][k];
          }
        }
    }
  }
}


//-----------------------------------------------------------------------------
// source term

void computeSourceTerm(const HeatEquation& physics, Real t, DiscVectorPtr rhs)
{
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
        rhs_arr[el][i] -= basis_vals.getValue(i, q) * weight * src_vals[q] / detJ[el][q];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Neumann BC

void computeNeumannBC(const HeatEquation& physics, const Real t, DiscVectorPtr u, DiscVectorPtr rhs)
{
  const auto& neumann_bcs = physics.getNeumannBCs();
  for (auto& bc : neumann_bcs)
  {
    std::cout << "doing Neumann BC " << bc->getSurfDisc()->getIdx() << std::endl;
    computeNeumannBC(bc, u, t, rhs);
  }
}

void computeNeumannBC(NeumannBCPtr bc, DiscVectorPtr u, const Real t, DiscVectorPtr rhs)
{
  // Note: we currently don't require that a given surface has a single volume as
  //       as its upward adjacency, so we have to look up a different volume disc
  //       for every face

  
  auto surf = bc->getSurfDisc();
  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::vector<Real> flux_vals(surf->getNumQuadPtsPerFace());
  Quadrature& quad = surf->quad;
  BasisVals2D basis(surf->face_group.getTPMapperSol(), quad.getPoints(), surf->face_group.getFaceNodesSol(), surf->face_group.ref_el_sol);

  for (int face=0; face < 6; ++face)
  {
    std::cout << "nodemap for face " << face << std::endl;
    for (int i=0; i < 4; ++i)
      std::cout << surf->face_group.getFaceNodesSol()[face][i] << ", ";

    std::cout << std::endl;
  }
  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    auto& face_spec = surf->face_group.faces[face];
    auto& u_arr = u->getArray(face_spec.vol_group);
    auto& res_arr = rhs->getArray(face_spec.vol_group);

    auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];
    auto res_el = res_arr[boost::indices[face_spec.el_group][range()]];

    surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);
    bc->getValue(face, t, u_quad.data(), flux_vals.data());
    ArrayType<Real, 1> face_sum(boost::extents[surf->getNumSolPtsPerFace()]);
    for (int k=0; k < surf->getNumSolPtsPerFace(); ++k)
      face_sum[k] = 0;


    for (int i=0; i < quad.getNumPoints(); ++i)
      for (int j=0; j < quad.getNumPoints(); ++j)
      {
        int node    = surf->quad_tp_nodemap[i][j];
        Real weight = quad.getWeight(i) * quad.getWeight(j);
        Real area   = std::sqrt(surf->normals[face][node][0] * surf->normals[face][node][0] +
                             surf->normals[face][node][1] * surf->normals[face][node][1] +
                             surf->normals[face][node][2] * surf->normals[face][node][2]);
        Real val = weight * area * flux_vals[node];
        //std::cout << "node   = " << node << ", flux val = " << flux_vals[node] << std::endl;
        //std::cout << "weight = " << weight << std::endl;
        //std::cout << "area   = " << area << std::endl;
        //std::cout << "val    = " << val << std::endl;

        for (int k=0; k < surf->getNumSolPtsPerFace(); ++k)
        {
          int node_sol = surf->face_group.getFaceNodesSol()[face_spec.face][k];
          //std::cout << "adding to el, node = " << face_spec.el_group << ", " << node_sol << std::endl;
          res_arr[face_spec.el_group][node_sol] -= basis.getValue(face_spec.face, k, i, j) * val;
          //std::cout << "contrib = " << basis.getValue(face_spec.face, k, i, j) * val << std::endl;
          face_sum[k] -= basis.getValue(face_spec.face, k, i, j) * val;
        }          
      }

    std::cout << "face (" << face_spec.el_group << ", " << face_spec.face << ") contribution sum: " << std::endl;
    for (int k=0; k < surf->getNumSolPtsPerFace(); ++k)
      std::cout << "node " << surf->face_group.getFaceNodesSol()[face_spec.face][k] << ": " << face_sum[k] << std::endl;
  }
}

void printArray(DiscVectorPtr vec)
{
  Real val_sum = 0;
  auto disc = vec->getDisc();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    std::cout << "Block " << i << std::endl;
    auto& arr_i = vec->getArray(i);

    for (unsigned int el=0; el < arr_i.shape()[0]; ++el)
    {
      Real el_sum = 0;
      for (unsigned int j=0; j < arr_i.shape()[1]; ++j)
      {
        std::cout << "arr(" << el << ", " << j << ") = " << arr_i[el][j] << std::endl;
        val_sum += arr_i[el][j];
        el_sum += arr_i[el][j];
      }

      std::cout << "el_sum = " << el_sum << std::endl;
    }
  }

  std::cout << "val_sum = " << val_sum << std::endl;
}


void printVector(DiscVectorPtr vec)
{
  std::cout << "printing vector form" << std::endl;
  vec->syncArrayToVector();

  auto disc = vec->getDisc();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    std::cout << "Block " << i << std::endl;
    auto vol_disc  = disc->getVolDisc(i);
    auto& dofs     = disc->getDofNumbering()->getDofs(i);
    auto& vec_vals = vec->getVector();
    std::cout << "dofs.shape = " << dofs.shape()[0] << ", " << dofs.shape()[1] << std::endl;
    std::cout << "vec_vals.shape = " << vec_vals.shape()[0] << std::endl;
    ArrayType<Real, 2> sol_coords(boost::extents[vol_disc->getNumSolPtsPerElement()][3]);

    for (int el=0; el < vol_disc->getNumElems(); ++el)
    {
      vol_disc->getVolumeSolCoords(el, sol_coords);
      for (unsigned int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
      {
        //std::cout << "el = " << el << ", j = " << j << std::endl;
        auto dof = dofs[el][j];
        if (disc->getDofNumbering()->isDofActive(dof))
        {
          //std::cout << "dof = " << dof << std::endl;
          std::cout << "arr(" << el << ", " << j << ") = " << vec_vals[dof] << std::endl;
          std::cout << "  coords = " << sol_coords[j][0] << ", " << sol_coords[j][1] << ", " << sol_coords[j][2] << std::endl;
        }
      }
    }
  }
}


}  // namespace 