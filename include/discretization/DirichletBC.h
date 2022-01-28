#ifndef DIRICHLETBC_H
#define DIRICHLETBC_H

#include "discretization/surface_discretization.h"
#include "discretization/boundary_condition.h"
#include "discretization/disc_vector.h"
#include "utils/error_handling.h"

class DirichletBC : public BoundaryCondition
{
  public:
    DirichletBC (SurfDiscPtr surf) : BoundaryCondition(surf)
    {
      assertAlways(surf->face_group.getIsDirichlet(), "Cannot define Dirichlet BC on non-Dirichlet surface");
    }

    virtual ~DirichletBC() {}

    // gets the prescribed value at the solution nodes on the face.
    // it would be better if the output argument were a boost array type,
    // but template functions cant be virtual
    virtual void getValue(const Index face, const Real t,  Real* vals) = 0;

};

using DirichletBCPtr = std::shared_ptr<DirichletBC>;

// sets the prescribed values in the array in disc_vec
void applyDirichletValues(DirichletBCPtr bc, const Real t, DiscVectorPtr disc_vec);

#endif
