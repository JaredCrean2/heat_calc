#ifndef TEST_MESH_HELPER_H
#define TEST_MESH_HELPER_H

#include <memory>
#include "mesh/mesh.h"
#include "mesh/mesh_input.h"
#include "discretization/surface_discretization.h"
#include "mesh/mesh_generator.h"
#include "quadrature.h"


Mesh::MeshSpec getStandardMeshSpec();

std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const Mesh::MeshSpec& spec = getStandardMeshSpec());


class StandardMeshSetup
{
  public:
    virtual ~StandardMeshSetup() {}

    // function must be called before first usage of any fields
    void setup(const int quad_degree = 5)
    {
      spec = getStandardMeshSpec();
      spec.nx = 4; spec.ny = 5, spec.nz = 6;

      //spec.nx = 1; spec.ny = 1, spec.nz = 1;
      //spec.xmin= 0; spec.ymin = 0; spec.zmin = 0;
      //spec.xmax = 1; spec.ymax = 1; spec.zmax = 1;
      setup(quad_degree, spec);
    }


    void setup(const int quad_degree, const Mesh::MeshSpec& spec)
    {
      quad = getGaussianQuadrature(quad_degree);

      mesh_dim_mins = {spec.xmin, spec.ymin, spec.zmin};
      mesh_dim_maxs = {spec.xmax, spec.ymax, spec.zmax};
      mesh = makeStandardMesh(spec);
 
      Mesh::VolumeGroup& vol_group = mesh->getElements(0);
      vol_disc = std::make_shared<VolumeDiscretization>(vol_group, quad);
      std::vector<std::shared_ptr<VolumeDiscretization>> vol_discs{vol_disc};
      
      for (Index i=0; i < mesh->getNumSurfaces(); ++i)
      {
        auto surf_i = std::make_shared<SurfaceDiscretization>(mesh->getFaces(i), quad, vol_discs); 
        surf_discs.push_back(surf_i);
      }
    }

    Mesh::MeshSpec spec;
    std::array<Real, 3> mesh_dim_mins;
    std::array<Real, 3> mesh_dim_maxs;
    std::shared_ptr<Mesh::MeshCG> mesh;
    Quadrature quad{getGaussianQuadrature(0)};
    VolDiscPtr vol_disc;
    std::vector<SurfDiscPtr> surf_discs;
};

#endif
