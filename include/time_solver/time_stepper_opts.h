#ifndef TIMESOLVERS_TIME_STEPPER_OPTS_H
#define TIMESOLVERS_TIME_STEPPER_OPTS_H

#include "ProjectDefs.h"
#include "timestep_controller.h"
#include "linear_system/large_matrix_factory.h"

namespace timesolvers {

struct TimeStepperOpts
{
  Real t_start = 0;
  Real t_end   = 0;
  std::shared_ptr<TimestepController> timestep_controller;
  linear_system::LargeMatrixType mat_type = linear_system::LargeMatrixType::Unknown;
  std::shared_ptr<linear_system::LargeMatrixOpts> matrix_opts = nullptr;
  Real nonlinear_abs_tol = -1;
  Real nonlinear_rel_tol = -1;
  int nonlinear_itermax  = -1;

  // Defines how the auxiliary equations are solved.  If false, a block system is
  // created that is iterated to convergence with Gauss-Seidel.  If true, a single
  // linear system is created with denser rows/columns appended to it for the
  // auxiliary equations
  bool solve_auxiliary_equations_combined_system = false;

  // if true, allocate 2 additional matrices, one for the mass matrix and
  // one for the linear (solution independent) part of the Jacobian.
  // This may speed the computation of the Jacobian at the expense of memory.
  bool precompute_linear_jacobian = false;
};

void checkTimeStepperOpts(const TimeStepperOpts& opts, bool check_implicit=true);


}

#endif