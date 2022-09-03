#ifndef HEAT_EQUATION_HELPER_FUNCS_H
#define HEAT_EQUATION_HELPER_FUNCS_H

#include "ProjectDefs.h"
#include "physics/heat/basis_vals.h"
#include <cassert>

namespace Heat {

template <typename T>
void zeroArray(ArrayType<T, 1>& arr)
{
  for (unsigned int i=0; i < arr.shape()[0]; ++i)
    arr[i] = 0;
}

template <typename T>
void zeroMatrix(ArrayType<T, 2>& arr)
{
  for (unsigned int i=0; i < arr.shape()[0]; ++i)
    for (unsigned j=0; j < arr.shape()[1]; ++j)
      arr[i][j] = 0;
}
/*
inline void computedNdx(const BasisVals& basis_vals, const ArrayType<Real, 4>& dxidx, int el, ArrayType<Real, 3>& dN_dx)
{
  auto numQuadPtsPerElement = dxidx.shape()[1];
  auto numSolPtsPerElement = dN_dx.shape()[0];
  assert(dN_dx.shape()[1] == numQuadPtsPerElement);
  assert(dN_dx.shape()[2] == 3);
  using Int = decltype(numQuadPtsPerElement);

  Real dNi_dxi[3];
  for (Int k=0; k < numQuadPtsPerElement; ++k)
    for (Int i=0; i < numSolPtsPerElement; ++i)
    {
      basis_vals.getDerivs(i, k, dNi_dxi);  // TODO: this is in the wrong order
      for (Int d1=0; d1 < 3; ++d1)
      {
        Real dN_dx_d1 = 0;
        for (Int d2=0; d2 < 3; ++d2)
          dN_dx_d1 += dxidx[el][k][d2][d1] * dNi_dxi[d2];  // TODO: store dxidx transpose to make this better
        dN_dx[i][k][d1] = dN_dx_d1;  // TODO: this is in the wrong order
      }
    }
}
*/

inline void computedNdx(const BasisVals& basis_vals, const ArrayType<Real, 4>& dxidx, int el, ArrayType<Real, 3>& dN_dx)
{
  auto numQuadPtsPerElement = dxidx.shape()[1];
  auto numSolPtsPerElement = dN_dx.shape()[0];
  assert(dN_dx.shape()[1] == numQuadPtsPerElement);
  assert(dN_dx.shape()[2] == 3);
  using Int = decltype(numQuadPtsPerElement);

  for (int i=0; i < numSolPtsPerElement; ++i)
    for (int k=0; k < numQuadPtsPerElement; ++k) 
      for (int d=0; d < 3; ++d)
        dN_dx[i][k][d] = 0;

  Real dNi_dxi[3];
  for (Int i=0; i < numSolPtsPerElement; ++i)
    #pragma clang loop vectorize(enable)
    for (Int k=0; k < numQuadPtsPerElement; ++k)
    {
      basis_vals.getDerivs(i, k, dNi_dxi);  // TODO: BasisVals uses a compressed representation, which makes
                                            //       this look like an indexed load
      #pragma clang loop vectorize(enable)
      for (Int d1=0; d1 < 3; ++d1)
      {
        //Real dN_dx_d1 = 0;
        #pragma clang loop vectorize(enable)
        for (Int d2=0; d2 < 3; ++d2)
          dN_dx[i][k][d1] += dxidx[el][k][d2][d1] * dNi_dxi[d2];
          //dN_dx_d1 += dxidx[el][k][d2][d1] * dNi_dxi[d2];  // TODO: store dxidx transpose to make this better
        //dN_dx[i][k][d1] = dN_dx_d1;  // TODO: this is in the wrong order
      }
    }
}

}

#endif 