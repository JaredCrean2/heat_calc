#include "physics/heat/volume_term.h"
#include "physics/heat/basis_vals.h"
#include "physics/heat/helper_funcs.h"

namespace Heat {

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
            rhs_arr[el][i] -= alpha * dNi_dx[d] * weight * dNj_dx[d] * u_arr[el][j] / detJ[el][k];
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
  Real alpha = params.kappa;
  auto& tp_mapper_sol = vol_disc->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_disc->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);
  ArrayType<Real, 3> dN_dx(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumQuadPtsPerElement()][3]);
  ArrayType<Real, 2> du_dx(boost::extents[vol_disc->getNumQuadPtsPerElement()][3]);

  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    std::cout << "\nel = " << el << std::endl;
    zeroMatrix(du_dx);
    computedNdx(basis_vals, dxidx, el, dN_dx);

    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      std::cout << "node " << i << ", u = " << u_arr[el][i] << std::endl;

    //TODO: is there a tensor-product form of this?
    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
        for (int d=0; d < 3; ++d)
          du_dx[k][d] += dN_dx[i][k][d] * u_arr[el][i];

    for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
      std::cout << "at quad pt " << k << ", du/dx = " << du_dx[k][0] << ", " << du_dx[k][1] << ", " << du_dx[k][2] << std::endl;

    for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
      {
        // compute dN_i/dx * weights * dN_j/dx * |J| * T_j
        auto& rev_nodemap = basis_vals.getRevNodemapOut();
        int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
        Real weight = vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);
        auto fac = alpha * weight / detJ[el][k];
        for (int d=0; d < 3; ++d)
          rhs_arr[el][i] -= fac * dN_dx[i][k][d] * du_dx[k][d];
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
    //computeVolumeTerm2Jac(vol_disc, params, u_arr, assembler);
    computeVolumeTerm3Jac(vol_disc, params, u_arr, assembler);

  }
}


void computeVolumeTerm2Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                           linear_system::AssemblerPtr assembler)
{

  auto& dxidx = vol_disc->dxidx;
  auto& detJ  = vol_disc->detJ;
  Real alpha = params.kappa;
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
          //TODO: store 1/detJ in an array to avoid the division here
          Real weight = alpha * vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k) / detJ[el][k];
          for (int d=0; d < 3; ++d)
            dR_du[i][j] -= dN_dx[i][k][d] * weight * dN_dx[j][k][d]; // / detJ[el][k];
        } 
      }

    assembler->assembleVolume(vol_disc->getIdx(), el, dR_du);
  }
}

void computeVolumeTerm3Jac_element(const int numSolPtsPerElement, const int numQuadPtsPerElement, 
                                          const Real alpha,
                                          const LocalIndex* const BOOST_RESTRICT rev_nodemap,
                                          const Real* const BOOST_RESTRICT weights,
                                          const Real* const BOOST_RESTRICT dN_dx,
                                          const Real* const BOOST_RESTRICT detJInv,
                                          Real* const BOOST_RESTRICT dR_du)
{

  #pragma clang loop vectorize(enable) interleave(enable)
  for (int i=0; i < numSolPtsPerElement; ++i)
    #pragma clang loop vectorize(enable) interleave(enable)
    for (int j=0; j < numSolPtsPerElement; ++j)
    {
      const int dR_du_idx = i * numSolPtsPerElement + j;
      //dR_du[dR_du_idx] = 0;
      Real val = 0;

      #pragma clang loop vectorize(enable) interleave(enable)
      for (int k=0; k < numQuadPtsPerElement; ++k)
      {
        int nodemap_idx = 3*k;
        int k_i = rev_nodemap[nodemap_idx]; int k_j = rev_nodemap[nodemap_idx + 1]; int k_k = rev_nodemap[nodemap_idx + 2];
        Real weight = alpha * weights[k_i] * weights[k_j] * weights[k_k] * detJInv[k];
        int dN_dx_idxi = i * numQuadPtsPerElement * 3 + k * 3;
        int dN_dx_idxj = j * numQuadPtsPerElement * 3 + k * 3;
        #pragma clang loop vectorize(enable) interleave(enable)
        for (int d=0; d < 3; ++d)
          val -= dN_dx[dN_dx_idxi + d] * weight * dN_dx[dN_dx_idxj + d]; // / detJ[el][k];
          //dR_du[dR_du_idx] -= dN_dx[dN_dx_idxi + d] * weight * dN_dx[dN_dx_idxj + d]; // / detJ[el][k];
      } 
      dR_du[dR_du_idx] = val;
    }
}

void computeVolumeTerm3Jac_element2(const int numSolPtsPerElement, const int numQuadPtsPerElement, 
                                          const Real alpha,
                                          const Real* const BOOST_RESTRICT weights,
                                          const Real* const BOOST_RESTRICT dN_dx,
                                          const Real* const BOOST_RESTRICT detJInv,
                                          Real* const BOOST_RESTRICT dR_du)
{

  #pragma clang loop vectorize(enable) interleave(enable)
  for (int i=0; i < numSolPtsPerElement; ++i)
    #pragma clang loop vectorize(enable) interleave(enable)
    for (int j=0; j < numSolPtsPerElement; ++j)
    {
      const int dR_du_idx = i * numSolPtsPerElement + j;
      //dR_du[dR_du_idx] = 0;
      Real val = 0;

      //#pragma clang loop vectorize(enable) interleave(enable)
      #pragma clang loop vectorize_width(4) interleave_count(4)
      for (int k=0; k < numQuadPtsPerElement; ++k)
      {
        Real weight = alpha * weights[k] * detJInv[k];
        int dN_dx_idxi = i * numQuadPtsPerElement * 3 + k * 3;
        int dN_dx_idxj = j * numQuadPtsPerElement * 3 + k * 3;
        #pragma clang loop vectorize(enable) interleave(enable)
        for (int d=0; d < 3; ++d)
          val -= dN_dx[dN_dx_idxi + d] * weight * dN_dx[dN_dx_idxj + d]; // / detJ[el][k];
          //dR_du[dR_du_idx] -= dN_dx[dN_dx_idxi + d] * weight * dN_dx[dN_dx_idxj + d]; // / detJ[el][k];
      } 
      dR_du[dR_du_idx] = val;
    }
}

void computeVolumeTerm3Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                           linear_system::AssemblerPtr assembler)
{

  auto& dxidx = vol_disc->dxidx;
  auto& detJInv  = vol_disc->detJInv;
  Real alpha = params.kappa;
  auto& tp_mapper_sol = vol_disc->vol_group.getTPMapperSol();
  Mesh::TensorProductMapper tp_mapper_quad(vol_disc->quad.getPoints());
  BasisVals basis_vals(tp_mapper_sol, tp_mapper_quad);
  const auto& rev_nodemap = basis_vals.getRevNodemapOut();

  ArrayType<Real, 3> dN_dx(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumQuadPtsPerElement()][3]);
  ArrayType<Real, 2> dR_du(boost::extents[vol_disc->getNumSolPtsPerElement()][vol_disc->getNumSolPtsPerElement()]);
  std::vector<Real> quad_weights(vol_disc->getNumQuadPtsPerElement());
  for (int k=0; k < vol_disc->getNumQuadPtsPerElement(); ++k)
  {
    int k_i = rev_nodemap[k][0]; int k_j = rev_nodemap[k][1]; int k_k = rev_nodemap[k][2];
    quad_weights[k] =  vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);
  }

  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    computedNdx(basis_vals, dxidx, el, dN_dx);
    //computeVolumeTerm3Jac_element(vol_disc->getNumSolPtsPerElement(), vol_disc->getNumQuadPtsPerElement(), alpha,
    //                              &(rev_nodemap[0][0]), vol_disc->quad.getWeights().data(), &(dN_dx[0][0][0]), &(detJInv[el][0]),
    //                              &(dR_du[0][0]));
    computeVolumeTerm3Jac_element2(vol_disc->getNumSolPtsPerElement(), vol_disc->getNumQuadPtsPerElement(), alpha,
                                    quad_weights.data(), &(dN_dx[0][0][0]), &(detJInv[el][0]),
                                    &(dR_du[0][0]));
    assembler->assembleVolume(vol_disc->getIdx(), el, dR_du);
  }
}


}
