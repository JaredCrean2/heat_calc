#ifndef NEUMANNBC_H
#define NEUMANNBC_H

#include "discretization/surface_discretization.h"
#include "discretization/boundary_condition.h"
#include "error_handling.h"

class NeumannBC : public BoundaryCondition
{
  public:
    NeumannBC (SurfDiscPtr surf, bool is_nonlinear) : 
      BoundaryCondition(surf),
      m_is_nonlinear(is_nonlinear)
    {
      assertAlways(!(surf->face_group.getIsDirichlet()), "Cannot define Neumann BC on a Dirichlet surface");
    }

    virtual ~NeumannBC() {}

    bool isNonlinear() const { return m_is_nonlinear; }

    // gets the prescribed flux value at the face quadrature nodes
    // it would be better if the output argument were a boost array type,
    // but template functions cant be virtual
    // flux_vals is 3 x numQuadPointsPerFace, giving the flux in each cartesian direction at each point
    virtual void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) = 0;

    // gets derivative of flux vals wrt sol_vals at each node.
    // flux_vals_deriv is 3 x numQuadPtsPerFace
    virtual void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv)
    {
      assert(!m_is_nonlinear);
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
        for (int d=0; d < 3; ++d)
          flux_vals_deriv[m_surf->getNumQuadPtsPerFace() * d + i] = 0;
    }

    std::array<Real, 3> getUnitNormal(int face, int quad_pt)
    {
      std::array<Real, 3> normal{m_surf->normals[face][quad_pt][0], m_surf->normals[face][quad_pt][1], m_surf->normals[face][quad_pt][2]};
      return normal / std::sqrt(dot(normal, normal));
    }

  private:
    bool m_is_nonlinear;

};

using NeumannBCPtr = std::shared_ptr<NeumannBC>;

#endif
