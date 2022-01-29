#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "discretization/surface_discretization.h"

class BoundaryCondition
{
  public:
    BoundaryCondition (SurfDiscPtr surf) : m_surf(surf) {}

    virtual ~BoundaryCondition() {}

    SurfDiscPtr getSurfDisc() const { return m_surf; }

  protected:
    
    // gets coordinates of all solution nodes
    void getSolNodeCoords(const Index face, ArrayType<Real, 2>& coords)
    {
      getSurfDisc()->getFaceSolCoords(face, coords);
    }

    // gets coordinates of all quadrature points
    void getQuadNodeCoords(const Index face, ArrayType<Real, 2>& coords)
    {
      getSurfDisc()->getFaceQuadCoords(face, coords);
    }

    SurfDiscPtr m_surf;
};

using BoundaryConditionPtr = std::shared_ptr<BoundaryCondition>;

#endif
