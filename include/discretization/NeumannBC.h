#ifndef DIRICHLETBC_H
#define DIRICHLETBC_H

#include "discretization/surface_discretization.h"

class NeumannBC
{
  public:
    NeumannBC (SurfDiscPtr surf) : m_surf(surf) {}

    virtual ~NeumannBC() {}

    // it would be better if the output argument were a boost array type,
    // but template functions cant be virtual
    virtual void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) const = 0;

    SurfDiscPtr getSurfDisc() const { return m_surf; }

  protected:

    void getNodeCoords(const Index face, const LocalIndex node,
                       Real coords[3])
    {
      auto& surf_group = getSurfDisc()->face_group;
      auto& face_spec = surf_group.faces[face];
      auto& vol_group = getSurfDisc()->volume_groups[face_spec.vol_group];
      auto& localidx = vol_group.nodemap_coord[face][node];

      for (int i=0; i < 3; ++i)
        coords[i] = vol_group.coords[face_spec.el_group][localidx][d];
    }

    SurfDiscPtr m_surf;
};

using NeumannBCPtr = std::shared_ptr<NeumannBC>;

#endif
