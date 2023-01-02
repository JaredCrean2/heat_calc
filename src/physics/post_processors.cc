#include "physics/post_processors.h"
#include "physics/heat/bc_defs.h"

namespace physics {

std::vector<double> PostProcessorBCFlux::getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
{
  if (!u->isArrayCurrent())
    u->syncArrayToVector();

  SurfDiscPtr surf = m_bc->getSurfDisc();
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


    m_bc->getValue(i, t, u_quad.data(), flux_vals.data());
    for (int j=0; j < surf->getNumQuadPtsPerFace(); ++j)
      for (int d=0; d < 3; ++d)
        flux_vals_array[j][d] = flux_vals[j + surf->getNumQuadPtsPerFace() * d];
      
    total_flux += integrateFaceVector(surf, i, flux_vals_array);
  }

  return {total_flux};
}

}