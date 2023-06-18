#include "physics/post_processors.h"
#include "ProjectDefs.h"
#include "physics/heat/bc_defs.h"

namespace physics {

Real integrateBoundaryFlux(NeumannBCPtr bc, DiscVectorPtr u, double t, MPI_Comm comm)
{
  if (!u->isArrayCurrent())
    u->syncVectorToArray();

  SurfDiscPtr surf = bc->getSurfDisc();
  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::vector<Real> flux_vals(3*surf->getNumQuadPtsPerFace());
  ArrayType<Real, 2> flux_vals_array(boost::extents[surf->getNumQuadPtsPerFace()][3]);

  Real total_flux = 0;
  for (int i=0; i < surf->getNumFaces(); ++i)
  {
    const Mesh::FaceSpec& spec = surf->face_group.faces[i];
    auto& u_arr = u->getArray(spec.vol_group);
    auto u_i = u_arr[boost::indices[spec.el_group][range()]];
    surf->interp_vsq_flat[spec.face].interpolateVals(u_i, u_quad);


    bc->getValue(i, t, u_quad.data(), flux_vals.data());
    for (int j=0; j < surf->getNumQuadPtsPerFace(); ++j)
      for (int d=0; d < 3; ++d)
        flux_vals_array[j][d] = flux_vals[j + surf->getNumQuadPtsPerFace() * d];
      
    total_flux += surf->getFaceWeight(i) * integrateFaceVector(surf, i, flux_vals_array);
  }

  Real val_global;
  MPI_Allreduce(&total_flux, &val_global, 1, REAL_MPI_DATATYPE, MPI_SUM, comm);

  return val_global;
}

//-----------------------------------------------------------------------------
// PostProcessorBCFlux

std::vector<double> PostProcessorBCFlux::getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
{
  return {integrateBoundaryFlux(m_bc, u, t, m_comm) };
}

//-----------------------------------------------------------------------------
// PostProcessorAirWindSkyBCFlux

std::vector<double> PostProcessorAirWindSkyBCFlux::getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
{
  Real interior_air_temp = u_aux->getVector(1)[0];
  m_heat_eqn_solar->setTimeParameters(t, interior_air_temp);
  return PostProcessorBCFlux::getValues(u, u_aux, t);
}

//-----------------------------------------------------------------------------
// PostProcessorCombinedAirWindSkyBCFlux

int PostProcessorCombinedAirWindSkyBCFlux::numValues() const { return m_bc->getNumBCs() + 1; }

std::vector<std::string> PostProcessorCombinedAirWindSkyBCFlux::getNames() const
{
  std::vector<std::string> names;
  for (int i=0; i < m_bc->getNumBCs(); ++i)
    names.push_back(m_name_prefix + "_" + m_bc->getBC(i)->getName());

  names.push_back(m_name_prefix);

  return names;
}

std::vector<double> PostProcessorCombinedAirWindSkyBCFlux::getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
{
  Real interior_air_temp = u_aux->getVector(1)[0];
  m_heat_eqn_solar->setTimeParameters(t, interior_air_temp);

  std::vector<double> vals;
  Real val_sum = 0.0;
  for (int i=0; i < m_bc->getNumBCs(); ++i)
  {
    vals.push_back(integrateBoundaryFlux(m_bc->getBC(i), u, t, m_comm));
    val_sum += vals.back();
  }

  vals.push_back(val_sum);

  return vals;
}

}