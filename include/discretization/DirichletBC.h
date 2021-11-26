#ifndef DIRICHLETBC_H
#define DIRICHLETBC_H

#include "discretization/surface_discretization.h"
#include "discretization/boundary_condition.h"

class DirichletBC : public BoundaryCondition
{
  public:
    DirichletBC (SurfDiscPtr surf) : BoundaryCondition(surf) {}

    virtual ~DirichletBC() {}

    // it would be better if the output argument were a boost array type,
    // but template functions cant be virtual
    virtual void getValue(const Index face, const Real t,  Real* vals) const = 0;

};

using DirichletBCPtr = std::shared_ptr<DirichletBC>;

#endif
