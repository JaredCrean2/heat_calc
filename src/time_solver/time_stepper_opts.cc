#include "time_solver/time_stepper_opts.h"

namespace timesolvers {

void checkTimeStepperOpts(const TimeStepperOpts& opts, bool check_implicit)
{
  assertAlways(opts.t_end > opts.t_start, "t_end must be > t_start");
  if (check_implicit)
  {
    assertAlways(opts.mat_type != linear_system::LargeMatrixType::Unknown, "Matrix type cannot be unknown");
    assertAlways(!!opts.matrix_opts, "matrix_opts must be defined");
    assertAlways(opts.nonlinear_abs_tol > 0 || opts.nonlinear_rel_tol > 0, "Either nonlinear_abs_tol or nonlinear_rel_tol must be > 0");
    assertAlways(opts.nonlinear_itermax > 0, "nonlinear_itermax must be > 0");
  }
}
}