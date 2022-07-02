#ifndef HEAT_SOURCE_TERM_H
#define HEAT_SOURCE_TERM_H

#include "physics/heat/HeatEquation.h"

namespace Heat {

void computeSourceTerm(const HeatEquation& physics, Real t, DiscVectorPtr rhs);

void computeSourceTerm(const VolDiscPtr vol_disc, SourceTermPtr src, Real t,
                       ArrayType<Real, 2>& rhs_arr);

}

#endif
