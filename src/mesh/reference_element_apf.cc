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
  int degree = ref_el->getNumNodes(1) + 1;
  std::string name = std::string("LagrangeHex") + std::to_string(degree);

  std::array<std::shared_ptr<EntityShape>, 4> entityshapes{std::make_shared<EntityShapeRefEl>(ref_el, 0, 0),
                                                                std::make_shared<EntityShapeRefEl>(ref_el, 1, 0),
                                                                std::make_shared<EntityShapeRefEl>(ref_el, 2, 0),
                                                                std::make_shared<EntityShapeRefEl>(ref_el, 3, 0)};

  return std::make_shared<FieldShapeRefEl>(name, ref_el, entityshapes);
}
}