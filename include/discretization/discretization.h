#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

#include <vector>
#include <memory>

#include "discretization/dof_numbering.h"
//#include "mesh/mesh.h"
#include "discretization/surface_discretization.h"
#include "discretization/volume_discretization.h"

namespace Mesh
{
  class MeshCG;
}


class Discretization
{
  public:

    Discretization(std::shared_ptr<Mesh::MeshCG>, const int volume_quad_accuracy,
                   const int surface_quad_accuracy);

    int getNumVolDiscs() const { return m_vol_discs.size(); }

    VolDiscPtr getVolDisc(const int i) { return m_vol_discs.at(i); }

    int getNumSurfDiscs() const { return m_surf_discs.size(); }

    int getNumBCSurfDiscs() const;

    SurfDiscPtr getSurfDisc(const int i) { return m_surf_discs.at(i); }

    SurfDiscPtr getBCSurfDisc(const int i)
    { 
      assertAlways(i < getNumBCSurfDiscs(), "invalid index for BC surface discretization");
      return m_surf_discs.at(i); 
    }

    std::shared_ptr<Mesh::MeshCG> getMesh() const { return m_mesh; }

    DofNumberingPtr getDofNumbering() const {return m_dof_numbering; }

  private:
    std::shared_ptr<Mesh::MeshCG> m_mesh;
    std::vector<VolDiscPtr> m_vol_discs;
    std::vector<SurfDiscPtr> m_surf_discs;
    DofNumberingPtr m_dof_numbering;
};

using DiscPtr = std::shared_ptr<Discretization>;

#endif
