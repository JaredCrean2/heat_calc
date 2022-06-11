#ifndef TIME_SOLVER_NEWTON_RESULT_H
#define TIME_SOLVER_NEWTON_RESULT_H

#include "ProjectDefs.h"
#include <ostream>

namespace timesolvers {

// describes the result of a Newton solve (converged, diverged etc)
class NewtonResult
{
  public:
    NewtonResult(Real norm0, Real norm, Real niters,
                 Real abstol, Real reltol, int itermax) :
      m_norm0(norm0),
      m_norm(norm),
      m_niters(niters),
      m_abstol(abstol),
      m_reltol(reltol),
      m_itermax(itermax)
    {}

    bool isAbstolSatisfied() const { return m_norm < m_abstol; }
    bool isReltolSatisfied() const { return m_norm/m_norm0 < m_reltol; }
    bool isItermaxSatisfied() const { return m_niters <= m_itermax; }  //TODO: need to make sure this works when m_niters = m_itermax
    bool isConverged() const { return (isAbstolSatisfied() || isReltolSatisfied()) && isItermaxSatisfied(); }

    Real getAbsNorm() const { return m_norm; }
    Real getRelNorm() const { return m_norm/m_norm0; }

    int getNIters() const { return m_niters; }
    Real getAbsTol() const { return m_abstol; }
    Real getRelTol() const { return m_reltol; }

  private:
    Real m_norm0;  // initial residual norm
    Real m_norm;   // final residual norm
    Real m_niters;

    Real m_abstol;
    Real m_reltol;
    int m_itermax;
};

std::ostream& operator<<(std::ostream& os, const NewtonResult& result);

} // namespace
#endif