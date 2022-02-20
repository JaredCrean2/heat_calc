#include "mesh/reference_element_apf.h"


namespace apf {

namespace Impl {

void fail(const std::string& why)
{
  apf::fail(why.c_str());
}
}

std::shared_ptr<FieldShapeRefEl> getHexFieldShape(std::shared_ptr<reference_element::ReferenceElement> ref_el)
{
  static std::map<std::shared_ptr<reference_element::ReferenceElement>, std::shared_ptr<FieldShapeRefEl>> fieldshapes;
  static std::vector<int> degree_counts;

  if (fieldshapes.count(ref_el) > 0)
    return fieldshapes[ref_el];

  int degree = ref_el->getDegree();

  int degree_count;
  if (degree_counts.size() <= degree)
  {
    degree_count = 0;
    degree_counts.resize(degree+1, 0);
    degree_counts[degree] += 1;
  } else
  {
    degree_count = degree_counts[degree];
    degree_counts[degree]++;
  }

  std::string name = std::string("LagrangeHex") + std::to_string(degree) + "_copy" + std::to_string(degree_count);

  std::array<std::shared_ptr<EntityShape>, 4> entityshapes{std::make_shared<EntityShapeRefEl>(ref_el, 0, 0),
                                                           std::make_shared<EntityShapeRefEl>(ref_el, 1, 0),
                                                           std::make_shared<EntityShapeRefEl>(ref_el, 2, 0),
                                                           std::make_shared<EntityShapeRefEl>(ref_el, 3, 0)};

  auto fieldshape = std::make_shared<FieldShapeRefEl>(name, ref_el, entityshapes);

  fieldshapes[ref_el] = fieldshape;

  return fieldshape;
}
}