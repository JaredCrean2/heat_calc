#ifndef TIMESOLVERS_TIMESTEP_CONTROLLER_PIECEWISE_H
#define TIMESOLVERS_TIMESTEP_CONTROLLER_PIECEWISE_H

#include "error_handling.h"
#include "time_solver/timestep_controller.h"

namespace timesolvers {

// determine the timestep by piecewise linear interplation
// of the given points
class TimestepControllerPiecewise : public TimestepController
{
  public:

    struct TimestepPoint
    {
      Real t;
      Real delta_t;
    };

    TimestepControllerPiecewise(const std::vector<TimestepPoint>& pts) :
      m_pts(pts)
    {
      assertAlways(pts.size() >= 2, "must have at least 2 points for TimestepControllerPiecewise");
      sortData();
      checkInputPts();
    }


    Real getNextTimestep(Real t) override
    {
      return interpolateDeltaT(t);
    }

    void recordLastIteration(Real physics_rhs) override {}

  private:
    void sortData()
    {
      auto cmp = [](const TimestepPoint& lhs, const TimestepPoint& rhs)
      {
        return lhs.t < rhs.t;
      };

      std::sort(m_pts.begin(), m_pts.end(), cmp);
    }

    void checkInputPts()
    {
      for (int i=1; i < int(m_pts.size()); ++i)
        assertAlways(m_pts[i].t > m_pts[i-1].t, "time points must be distinct");
    }

    int getStartIdx(Real t)
    {
      if (t < m_pts[0].t)
        throw std::runtime_error("t < first time point");

      if (t > m_pts.back().t)
        throw std::runtime_error("t > last time point");

      for (int i=0; i < int(m_pts.size() - 1); ++i)
        if (t >= m_pts[i].t && t <= m_pts[i+1].t)
          return i;

      throw std::runtime_error("unable to determine time index");
    }

    Real interpolateDeltaT(Real t)
    {
      int idx = getStartIdx(t);

      TimestepPoint p1 = m_pts[idx];
      TimestepPoint p2 = m_pts[idx+1];

      Real fac = (t - p1.t)/(p2.t - p1.t);

      return p1.delta_t + fac * (p2.delta_t - p1.delta_t);
    }


    std::vector<TimestepPoint> m_pts;
};
}

#endif