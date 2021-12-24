#include "physics/PhysicsModel.h"
#include <stdexcept>

void PhysicsModel::checkInitialization()
{
  if (m_source_terms.size() != static_cast<size_t>(getDiscretization()->getNumVolDiscs()))
    throw std::runtime_error("number of source terms must equal number of volume discretizations");

}
