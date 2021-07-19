#ifndef DIRICHLETBC_H
#define DIRICHLETBC_H

#include "discretization/surface_discretization.h"

class DirichletBC
{
  public:
    DirichletBC (SurfDiscPtr surf) : m_surf(surf) {}

    virtual ~DirichletBC() {}

    // it would be better if the output argument were a boost array type,
    // but template functions cant be virtual
    virtual void getValue(const Index face, const Real t,  Real* vals) const = 0;

    SurfDiscPtr getSurfDisc() const { return m_surf; }

  protected:

    // gets coordinates of all nodes on a given face
    void getNodeCoords(const Index face, ArrayType<Real, 2>& coords)
    {
      auto& surf_group = getSurfDisc()->face_group;
      auto& face_spec = surf_group.faces[face];
      auto& vol_disc = getSurfDisc()->volume_discs[face_spec.vol_group];
      //auto& localidx = vol_disc.nodemap_coord[face][node];
      auto& mapper = vol_disc->vol_group.getTPMapperCoord();
      ArrayType<Real, 3> coords_tp(mapper.getTPShape());  //TODO: store this
      for (int d=0; d < 3; ++d)
      {
        auto coords_d = vol_disc->vol_group.coords[boost::indices[face_spec.el_group][range()][d]];
        auto coords_out = coords[boost::indices[range()][d]];
        mapper.mapToTP(coords_d, coords_tp);
        getSurfDisc()->interp_vcq_tp[face].interpolateVals(coords_tp, coords_out);
      }
    }

    SurfDiscPtr m_surf;
};

using DirichletBCPtr = std::shared_ptr<DirichletBC>;

#endif
