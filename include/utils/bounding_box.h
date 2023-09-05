#ifndef HEAT_CALC_UTILS_BOUNDING_BOX_H
#define HEAT_CALC_UTILS_BOUNDING_BOX_H

#include <array>
#include "ProjectDefs.h"
#include "error_handling.h"
#include <iostream>

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

    BoundingBox() :
      m_lower_left({MIN_REAL, MIN_REAL, MIN_REAL}),
      m_upper_right({MAX_REAL, MAX_REAL, MAX_REAL})
    {}

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

  friend std::ostream& operator<<(std::ostream& os, const BoundingBox& box);

};

inline std::ostream& operator<<(std::ostream& os, const BoundingBox& box)
{
  std::cout << "bounding box from " << box.m_lower_left[0] << ", " << box.m_lower_left[1] << ", " << box.m_lower_left[2] 
            << " to " << box.m_upper_right[0] << ", " << box.m_upper_right[1] << ", " << box.m_upper_right[2];
  return os;
}

}

#endif