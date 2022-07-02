#ifndef HEAT_NEUMANN_BC_H
#define HEAT_NEUMANN_BC_H

#include "physics/heat/HeatEquation.h"

namespace Heat {

void computeNeumannBC(const HeatEquation& physics, const Real t, DiscVectorPtr u, DiscVectorPtr rhs);

void computeNeumannBC(NeumannBCPtr bc, DiscVectorPtr u, const Real t, DiscVectorPtr rhs);

void computeNeumannBCJacobian(const HeatEquation& physics, DiscVectorPtr u, Real t, linear_system::AssemblerPtr assembler);

void computeNeumannBCJacobian(NeumannBCPtr bc, DiscVectorPtr u, Real t, linear_system::AssemblerPtr assembler);

}

#endif
