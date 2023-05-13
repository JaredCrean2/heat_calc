#ifndef DIRICHLETBC_H
#define DIRICHLETBC_H

#include "discretization/surface_discretization.h"
#include "discretization/boundary_condition.h"
#include "discretization/disc_vector.h"
#include "utils/error_handling.h"

class DirichletBC : public BoundaryCondition
{
  public:
    DirichletBC (SurfDiscPtr surf, bool is_unsteady) : 
      BoundaryCondition(surf),
       m_is_unsteady(is_unsteady)
    {
      assertAlways(surf->face_group.getIsDirichlet(), "Cannot define Dirichlet BC on non-Dirichlet surface");
    }

    virtual ~DirichletBC() {}

    bool getIsUnsteady() const { return m_is_unsteady; }

    // gets the prescribed value at the solution nodes on the face.
    // it would be better if the output argument were a boost array type,
    // but template functions cant be virtual
    virtual void getValue(const Index face, const Real t,  Real* vals) = 0;

    // returns the time derivative of getValue.  This function will only
    // be called if getIsUnsteady() is true
    virtual void getValueDt(const Index face, const Real t, Real* vals) = 0;

  private:
    bool m_is_unsteady;

};

using DirichletBCPtr = std::shared_ptr<DirichletBC>;

// sets the prescribed values in the array in disc_vec
void applyDirichletValues(DirichletBCPtr bc, const Real t, DiscVectorPtr disc_vec);

// copies values from the dirichlet nodes on the the block that has the dirichlet surface
// to the same nodes on other blocks
void updateDependentDirichletValues(DiscVectorPtr disc_vec);


namespace impl {

Real errorFunc(Real x, Real y, Real z, Real t);

Real zeroFunc(Real x, Real y, Real z, Real t);

using ErrorFuncType = Real (*)(Real, Real, Real, Real);
using ZeroFuncType  = Real (*)(Real, Real, Real, Real);


}

#endif
