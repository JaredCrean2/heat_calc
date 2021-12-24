#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "discretization/surface_discretization.h"

class BoundaryCondition
{
  public:
    BoundaryCondition (SurfDiscPtr surf) : m_surf(surf) {}

    virtual ~BoundaryCondition() {}

    SurfDiscPtr getSurfDisc() const { return m_surf; }

  protected:
    
    // gets coordinates of all solution nodes
    void getSolNodeCoords(const Index face, ArrayType<Real, 2>& coords)
    {
      auto& surf_group = getSurfDisc()->face_group;
      auto& face_spec = surf_group.faces[face];
      auto& vol_disc = getSurfDisc()->volume_discs[face_spec.vol_group];
      //auto& localidx = vol_disc.nodemap_coord[face][node];
      for (int d=0; d < 3; ++d)
      {
        auto coords_d = vol_disc->vol_group.coords[boost::indices[face_spec.el_group][range()][d]];
        auto coords_out = coords[boost::indices[range()][d]];
        getSurfDisc()->interp_vcs_flat[face_spec.face].interpolateVals(coords_d, coords_out);
      }
    }

    // gets coordinates of all quadrature points
    void getQuadNodeCoords(const Index face, ArrayType<Real, 2>& coords)
    {
      auto& surf_group = getSurfDisc()->face_group;
      auto& face_spec = surf_group.faces[face];
      auto& vol_disc = getSurfDisc()->volume_discs[face_spec.vol_group];
      //auto& localidx = vol_disc.nodemap_coord[face][node];
      for (int d=0; d < 3; ++d)
      {
        auto coords_d = vol_disc->vol_group.coords[boost::indices[face_spec.el_group][range()][d]];
        auto coords_out = coords[boost::indices[range()][d]];
        getSurfDisc()->interp_vcq_flat[face].interpolateVals(coords_d, coords_out);
      }
    }

    SurfDiscPtr m_surf;
};

using BoundaryConditionPtr = std::shared_ptr<BoundaryCondition>;

#endif
