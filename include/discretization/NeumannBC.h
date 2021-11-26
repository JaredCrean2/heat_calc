#ifndef NEUMANNBC_H
#define NEUMANNBC_H

#include "discretization/surface_discretization.h"
#include "discretization/boundary_condition.h"

class NeumannBC : public BoundaryCondition
{
  public:
    NeumannBC (SurfDiscPtr surf) : BoundaryCondition(surf) {}

    virtual ~NeumannBC() {}

    // it would be better if the output argument were a boost array type,
    // but template functions cant be virtual
    virtual void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) const = 0;

};

using NeumannBCPtr = std::shared_ptr<NeumannBC>;

#endif
