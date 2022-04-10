#include "physics/heat/HeatEquation.h"
#include "discretization/disc_vector.h"
#include "discretization/dof_numbering.h"
#include "physics/heat/basis_vals.h"
#include "physics/heat/helper_funcs.h"
#include "mesh/mesh.h"

#include "linear_system/assembler.h"

namespace Heat {


void HeatEquation::computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs)
{
  rhs->set(0);
  if (!(u->isArrayCurrent()))
  {
    //std::cout << "syncing u vector to array" << std::endl;
    u->syncVectorToArray();
  }

  //std::cout << "u array = " << std::endl;
  //printArray(u);

  //std::cout << "before setting dirichlet values, u = " << std::endl;
  //printArray(u);

  // apply Dirichlet values to array
  applyDirichletValues(*this, t, u);
  //std::cout << "after setting dirichlet values, u = " << std::endl;
  //printArray(u);

  // compute volume terms
  computeVolumeTerm(*this, u, rhs);

  //std::cout << "\nafter volume term " << std::endl;
  //printArray(rhs);

  // compute Neumann BC terms
  computeNeumannBC(*this, t, u, rhs);

  //std::cout << "\nafter Neumann term " << std::endl;
  //printArray(rhs);

  // compute source term
  computeSourceTerm(*this, t, rhs);

  //std::cout << "\nafter source term " << std::endl;
  //printArray(rhs);
  //printVector(rhs);
}


void HeatEquation::computeJacobian(DiscVectorPtr u, const Real t, linear_system::AssemblerPtr assembler)
{
  if (!u->isArrayCurrent())
    u->syncArrayToVector();

  applyDirichletValues(*this, t, u);

  computeVolumeJacobian(*this, u, assembler);

  // typical Neumann and source terms don't contribute to the Jacobian
}

void HeatEquation::applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out)
{
  if (!vec_in->isArrayCurrent())
    vec_in->syncVectorToArray();

  ::Heat::applyMassMatrix(*this, vec_in, vec_out);
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
// ApplyMassMatrix

void applyMassMatrix(const HeatEquation& physics, DiscVectorPtr vec_in, DiscVectorPtr vec_out)
{
  if (!vec_in->isArrayCurrent())
    vec_in->syncVectorToArray();

  auto disc = physics.getDiscretization();

  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto vol_disc = disc->getVolDisc(i);
    auto dofs     = disc->getDofNumbering();
    auto& arr_in  = vec_in->getArray(i);
    auto& arr_out = vec_out->getArray(i);

    applyMassMatrix(vol_disc, dofs, arr_in, arr_out);
  }

  vec_out->markArrayModified();
}

void applyMassMatrix(const VolDiscPtr vol_disc, const DofNumberingPtr dof_numbering,
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
      u_quad[k] = u_quad[k] * weight / detJ[el][k];
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

    //computeVolumeTerm(vol_disc, params, u_arr, rhs_arr);
    computeVolumeTerm2(vol_disc, params, u_arr, rhs_arr);
  }

  rhs->markArrayModified();
}

//TODO: is this still used?
std::array<Real, 3> interpolateXi(const BasisVals& basis_vals, int k)
{
  ArrayType<Real, 2> xi_vals(boost::extents[8][3]);
  for (int i=0; i < 8; ++i)
    for (int j=0; j < 3; ++j)
      xi_vals[i][j] = 0;

  
  xi_vals[1][0] = 1;  xi_vals[2][1] = 1;  xi_vals[4][2] = 1;
  xi_vals[2][0] = 1;  xi_vals[3][1] = 1;  xi_vals[5][2] = 1;
  xi_vals[5][0] = 1;  xi_vals[6][1] = 1;  xi_vals[6][2] = 1;
  xi_vals[6][0] = 1;  xi_vals[7][1] = 1;  xi_vals[7][2] = 1;

  std::array<Real, 3> xi_interp{0, 0, 0};
  for (int dim=0; dim < 3; ++dim)
    for (int j=0; j < 8; ++j)
    {
      xi_interp[dim] += basis_vals.getValue(j, k) * xi_vals[j][dim];
    }

  return xi_interp;
}


void computeVolumeTerm(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                       const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr)
{
  auto& dxidx = vol_disc->dxidx;
  auto& detJ  = vol_disc->detJ;
  Real alpha = params.kappa / (params.rho * params.Cp);
  auto& tp_mapper_sol = vol_disc->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_disc->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);
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


void computeVolumeTerm2(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                        const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr)
{
  auto& dxidx = vol_disc->dxidx;
  auto& detJ  = vol_disc->detJ;
  Real alpha = params.kappa / (params.rho * params.Cp);
  auto& tp_mapper_sol = vol_disc->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_disc->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);
  ArrayType<Real, 3> dN_dx(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumQuadPtsPerElement()][3]);
  ArrayType<Real, 2> du_dx(boost::extents[vol_disc->getNumQuadPtsPerElement()][3]);

  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    zeroMatrix(du_dx);
    computedNdx(basis_vals, dxidx, el, dN_dx);

    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
        for (int d=0; d < 3; ++d)
          du_dx[k][d] += dN_dx[i][k][d] * u_arr[el][i];

    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
      {
        // compute dN_i/dx * weights * dN_j/dx * |J| * T_j
        auto& rev_nodemap = basis_vals.getRevNodemapOut();
        int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
        Real weight = vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);
        for (int d=0; d < 3; ++d)
        {
          rhs_arr[el][i] += alpha * dN_dx[i][k][d] * weight * du_dx[k][d] / detJ[el][k];
        }
      }
  }
}

//-----------------------------------------------------------------------------
// volume term Jacobian

void computeVolumeJacobian(const HeatEquation& physics, DiscVectorPtr u, linear_system::AssemblerPtr assembler)
{
  auto disc = physics.getDiscretization();

  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto vol_disc = disc->getVolDisc(i);
    auto& params  = physics.getVolumeGroupParams(i);
    auto& u_arr   = u->getArray(i);

    //computeVolumeTerm(vol_disc, params, u_arr, rhs_arr);
    computeVolumeTerm2Jac(vol_disc, params, u_arr, assembler);
  }
}


void computeVolumeTerm2Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                            linear_system::AssemblerPtr assembler)
{

  auto& dxidx = vol_disc->dxidx;
  auto& detJ  = vol_disc->detJ;
  Real alpha = params.kappa / (params.rho * params.Cp);
  auto& tp_mapper_sol = vol_disc->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_disc->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);
  const auto& rev_nodemap = basis_vals.getRevNodemapOut();

  ArrayType<Real, 3> dN_dx(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumQuadPtsPerElement()][3]);
  ArrayType<Real, 2> dR_du(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumSolPtsPerElement()]);

  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    computedNdx(basis_vals, dxidx, el, dN_dx);

    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      for (int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
      {
        dR_du[i][j] = 0;

        for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
        {
          int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
          Real weight = vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);
          for (int d=0; d < 3; ++d)
            dR_du[i][j] += alpha * dN_dx[i][k][d] * weight * dN_dx[j][k][d] / detJ[el][k];
        } 
      }

    assembler->assembleVolume(vol_disc->getIdx(), el, dR_du);
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

//TODO: I don't think the density and specific heat capacity are correctly accounted for here

void computeNeumannBC(const HeatEquation& physics, const Real t, DiscVectorPtr u, DiscVectorPtr rhs)
{
  const auto& neumann_bcs = physics.getNeumannBCs();
  for (auto& bc : neumann_bcs)
    computeNeumannBC(bc, u, t, rhs);

  rhs->markArrayModified();
}

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
    auto& face_spec = surf->face_group.faces[face];
    auto& u_arr = u->getArray(face_spec.vol_group);
    auto& res_arr = rhs->getArray(face_spec.vol_group);

    auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];
    auto res_el = res_arr[boost::indices[face_spec.el_group][range()]];

    surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);
    bc->getValue(face, t, u_quad.data(), flux_vals.data());

    for (int i=0; i < quad.getNumPoints(); ++i)
      for (int j=0; j < quad.getNumPoints(); ++j)
      {
        int node    = surf->quad_tp_nodemap[i][j];
        Real weight = quad.getWeight(i) * quad.getWeight(j);
        Real flux_normal = 0;
        for (int d=0; d < 3; ++d)
          flux_normal += surf->normals[face][node][d] * flux_vals[node + surf->getNumQuadPtsPerFace() * d];

        Real val = weight * flux_normal;


        for (int k=0; k < surf->getNumSolPtsPerFace(); ++k)
        {
          int node_sol = surf->face_group.getFaceNodesSol()[face_spec.face][k];
          res_arr[face_spec.el_group][node_sol] -= basis.getValue(face_spec.face, k, i, j) * val;
        }          
      }
  }
}

void printArray(DiscVectorPtr vec)
{
  Real val_sum = 0;
  auto disc = vec->getDisc();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    std::cout << "Block " << i << std::endl;
    auto vol_disc  = disc->getVolDisc(i);
    ArrayType<Real, 2> sol_coords(boost::extents[vol_disc->getNumSolPtsPerElement()][3]);

    auto& arr_i = vec->getArray(i);

    for (unsigned int el=0; el < arr_i.shape()[0]; ++el)
    {
      vol_disc->getVolumeSolCoords(el, sol_coords);
      Real el_sum = 0;
      for (unsigned int j=0; j < arr_i.shape()[1]; ++j)
      {
        //std::cout << "arr(" << el << ", " << j << ") = " << arr_i[el][j] << std::endl;

        std::string coord_str = std::string("coords = ") + std::to_string(sol_coords[j][0]) + ", " + 
                                std::to_string(sol_coords[j][1]) + ", " + std::to_string(sol_coords[j][2]) + ", ";
        int nspaces = std::max(40 - coord_str.size(), size_t(1));
        std::string space_str(nspaces, ' ');
        std::cout << "el " << el << " " << coord_str << space_str << "arr(" << el << ", " << j << ") = " << arr_i[el][j] << std::endl;

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
          //std::cout << "arr(" << el << ", " << j << ") = " << vec_vals[dof] << std::endl;
          std::string coord_str = std::string("coords = ") + std::to_string(sol_coords[j][0]) + ", " + 
                                  std::to_string(sol_coords[j][1]) + ", " + std::to_string(sol_coords[j][2]) + ", ";
          int nspaces = std::max(40 - coord_str.size(), size_t(1));
          std::string space_str(nspaces, ' ');
          std::cout << "el " << el << " " << coord_str << space_str << "arr(" << el << ", " << j << ") = " << vec_vals[dof] << std::endl;

         // std::cout << "coords = " << sol_coords[j][0] << ", " << sol_coords[j][1] << ", " << sol_coords[j][2]
         //           << ",  arr(" << el << ", " << j << ") = " << vec_vals[dof] << std::endl;
        }
      }
    }
  }
}


}  // namespace 