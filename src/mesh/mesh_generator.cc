#include "mesh/mesh_generator.h"

namespace Mesh
{

Point identity(const Point& point)
{
  return {point.x, point.y, point.z};
}

} // namespace
