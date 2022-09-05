#ifndef TIME_SOLVER_NEWTON_RESULT_H
#define TIME_SOLVER_NEWTON_RESULT_H

#include "ProjectDefs.h"
#include <ostream>

namespace timesolvers {

//TODO: move to own file
struct NewtonOpts
{
  Real nonlinear_abs_tol = 1e-8;
  Real nonlinear_rel_tol = 1e-8;
  int nonlinear_itermax  = 10;

  // parameters for Gauss-Seidel iteration
  int linear_itermax = 1;
  Real linear_delta_u_tol = 1e-2;  // terminate Gauss-Seidel iteration if
                                   // the current iteration's contribution to
                                   // delta_u is less than the total delta_u
};


// describes the result of a Newton solve (converged, diverged etc)
class NewtonResult
{
  public:
    NewtonResult(Real norm0, Real norm, Real niters,
                 NewtonOpts opts) :
      m_norm0(norm0),
      m_norm(norm),
      m_niters(niters),
      m_opts(opts)
    {}

    bool isAbstolSatisfied() const { return m_norm < m_opts.nonlinear_abs_tol; }
    bool isReltolSatisfied() const { return m_norm/m_norm0 < m_opts.nonlinear_rel_tol; }
    bool isItermaxSatisfied() const { return m_niters <= m_opts.nonlinear_itermax; }  //TODO: need to make sure this works when m_niters = m_itermax
    bool isConverged() const { return (isAbstolSatisfied() || isReltolSatisfied()) && isItermaxSatisfied(); }

    Real getAbsNorm() const { return m_norm; }
    Real getRelNorm() const { return m_norm/m_norm0; }

    int getNIters() const { return m_niters; }
    Real getAbsTol() const { return m_opts.nonlinear_abs_tol; }
    Real getRelTol() const { return m_opts.nonlinear_rel_tol; }

  private:
    Real m_norm0;  // initial residual norm
    Real m_norm;   // final residual norm
    Real m_niters;

    NewtonOpts m_opts;
};

std::ostream& operator<<(std::ostream& os, const NewtonResult& result);

} // namespace
#endif