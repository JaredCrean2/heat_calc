#include "physics/heat/post_processor_interior.h"
#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/AuxiliaryEquationsSolar.h"


namespace Heat {

std::vector<Real> PostProcessorInterior::getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
{
  Real hvac_flux = m_air_temp->getHVACFlux();
  Real t_air = u_aux->getVector(1)[0];

  return {t_air, hvac_flux};
}

}  // namespace