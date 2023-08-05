#ifndef HEAT_CALC_UTILS_BOUNDING_BOX_H
#define HEAT_CALC_UTILS_BOUNDING_BOX_H

#include <array>
#include "ProjectDefs.h"
#include "error_handling.h"

namespace utils {

class BoundingBox
{
  public:
    BoundingBox(const std::array<Real, 3>& lower_left,
                const std::array<Real, 3>& upper_right) :
      m_lower_left(lower_left),
      m_upper_right(upper_right)
    {
      for (int i=0; i < 3; ++i)
        assertAlways(upper_right[i] > lower_left[i], "upper_right must be greater than lower_left");
    }

    bool contains(const std::array<Real, 3>& pt)
    {
      bool found = true;
      for (int i=0; i < 3; ++i)
        found = found && (pt[i] >= m_lower_left[i] && pt[i] <= m_upper_right[i]);

      return found;
    }

  private:
    std::array<Real, 3> m_lower_left;
    std::array<Real, 3> m_upper_right;
};

}

#endif