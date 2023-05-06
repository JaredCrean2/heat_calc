#ifndef HEAT_VOLUME_TERM_H
#define HEAT_VOLUME_TERM_H

#include "HeatEquation.h"

namespace Heat {

void computeVolumeTerm(const HeatEquation& physics, DiscVectorPtr u, DiscVectorPtr rhs);

void computeVolumeTerm(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                       const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);

void computeVolumeTerm2(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                        const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);

void computeVolumeJacobian(const HeatEquation& physics, DiscVectorPtr u, linear_system::AssemblerPtr assembler);

void computeVolumeTerm2Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                            linear_system::AssemblerPtr assembler);

void computeVolumeTerm3Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                           linear_system::AssemblerPtr assembler);

void computeVolumeTerm(const HeatEquation& physics, DiscVectorPtr u, DiscVectorPtr rhs);

void computeVolumeTerm(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                       const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);

void computeVolumeTerm2(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                        const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);

void computeVolumeTerm3(const VolDiscPtr vol_disc, const VolumeGroupParams& params,
                        const ArrayType<Real, 2>& u_arr, ArrayType<Real, 2>& rhs_arr);                        

void computeVolumeJacobian(const HeatEquation& physics, DiscVectorPtr u, linear_system::AssemblerPtr assembler);

void computeVolumeTerm2Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                            linear_system::AssemblerPtr assembler);

void computeVolumeTerm3Jac(const VolDiscPtr vol_disc, const VolumeGroupParams& params, const ArrayType<Real, 2> u_arr,
                           linear_system::AssemblerPtr assembler);

}

#endif
