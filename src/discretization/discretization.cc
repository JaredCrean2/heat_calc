#include "discretization/discretization.h"
#include "mesh/mesh.h"

Discretization::Discretization(std::shared_ptr<Mesh::MeshCG> mesh,
          const int volume_quad_accuracy, const int surface_quad_accuracy) :
  m_mesh(mesh),
  m_vol_discs(mesh->getNumVolumeGroups()),
  m_surf_discs(mesh->getNumSurfaces()),
  m_dof_numbering(std::make_shared<DofNumbering>(mesh))
{
  auto volume_quad = getGaussianQuadrature(volume_quad_accuracy);
  for (int i=0; i < mesh->getNumVolumeGroups(); ++i)
    m_vol_discs[i] = std::make_shared<VolumeDiscretization>(mesh->getElements(i), volume_quad);

  auto surf_quad = getGaussianQuadrature(surface_quad_accuracy);
  for (int i=0; i < mesh->getNumSurfaces(); ++i)
    m_surf_discs[i] = std::make_shared<SurfaceDiscretization>(mesh->getFaces(i), surf_quad, m_vol_discs);
}

int Discretization::getNumBCSurfDiscs() const { return m_mesh->getNumBCSurfaces(); }