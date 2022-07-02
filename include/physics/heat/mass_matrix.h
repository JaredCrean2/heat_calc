#ifndef HEAT_EQUATION_MASS_MATRIX
#define HEAT_EQUATION_MASS_MATRIX

#include "HeatEquation.h"

namespace Heat {

void applyMassMatrix(const HeatEquation& physics, DiscVectorPtr vec_in, DiscVectorPtr vec_out);

void applyMassMatrix(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const DofNumberingPtr dof_numbering,
                     const ArrayType<Real, 2>& arr_in, ArrayType<Real, 2>& arr_out);

void computeMassMatrix(const HeatEquation& physics, linear_system::AssemblerPtr assembler);

void computeMassMatrix(const VolDiscPtr vol_disc, const VolumeGroupParams& params, linear_system::AssemblerPtr assembler);

}
#endif
