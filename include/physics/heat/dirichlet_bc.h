#ifndef HEAT_DIRICHLET_BC_H
#define HEAT_DIRICHLET_BC_H

#include "physics/heat/HeatEquation.h"

namespace Heat {

void applyDirichletValues(const HeatEquation& physics, const Real t, DiscVectorPtr u);

void computeUnsteadyDirichletBC(const HeatEquation& physics, const Real t, DiscVectorPtr rhs);

void computeUnsteadyDirichletBC(DirichletBCPtr bc, const std::vector<VolDiscPtr>& vol_discs,
                                const std::vector<const VolumeGroupParams*>& vol_group_params,
                                const Real t, DiscVectorPtr rhs);

}

#endif
