#ifndef TEST_MESH_HELPER_H
#define TEST_MESH_HELPER_H

#include <memory>
#include "mesh/mesh.h"
#include "mesh/mesh_input.h"
#include "discretization/surface_discretization.h"
#include "discretization/discretization.h"
#include "mesh/mesh_generator.h"
#include "quadrature.h"


Mesh::MeshSpec getStandardMeshSpec();

std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const Mesh::MeshSpec& spec = getStandardMeshSpec(), int sol_degree=1,
                                const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true});

// creates a single block mesh with all dirichlet BCs
class StandardMeshBase
{
  public:
    virtual ~StandardMeshBase() {}

    // function must be called before first usage of any fields
    virtual void setup(const int quad_degree=5, int sol_degree=1, 
              const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true})
    {
      spec = getStandardMeshSpec();
      spec.nx = 4; spec.ny = 5, spec.nz = 6;

      //spec.nx = 2; spec.ny = 1, spec.nz = 1;
      //spec.xmin = 0; spec.ymin = 0; spec.zmin = 0;
      //spec.xmax = 1; spec.ymax = 1; spec.zmax = 1;
      setup(quad_degree, sol_degree, spec, is_surf_dirichlet);
    }


    virtual void setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& meshspec,
                       const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true})
    {
      spec = meshspec;
      quad = getGaussianQuadrature(quad_degree);

      mesh_dim_mins = {spec.xmin, spec.ymin, spec.zmin};
      mesh_dim_maxs = {spec.xmax, spec.ymax, spec.zmax};
      mesh = makeStandardMesh(spec, sol_degree, is_surf_dirichlet);
      

    }

    Mesh::MeshSpec spec;
    std::array<Real, 3> mesh_dim_mins;
    std::array<Real, 3> mesh_dim_maxs;
    std::shared_ptr<Mesh::MeshCG> mesh;
    Quadrature quad{getGaussianQuadrature(0)};
};


// creates mesh + volume disc and 1 surface disc for each face of the cube
class StandardMeshSetup : public StandardMeshBase
{
  public:
    virtual ~StandardMeshSetup() {}

    using StandardMeshBase::setup;

    virtual void setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& spec,
                       const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true}) override
    {
      StandardMeshBase::setup(quad_degree, sol_degree, spec, is_surf_dirichlet);
 
      Mesh::VolumeGroup& vol_group = mesh->getElements(0);
      vol_disc = std::make_shared<VolumeDiscretization>(vol_group, quad);
      std::vector<std::shared_ptr<VolumeDiscretization>> vol_discs{vol_disc};
      
      surf_discs.resize(0);
      for (Index i=0; i < mesh->getNumSurfaces(); ++i)
      {
        auto surf_i = std::make_shared<SurfaceDiscretization>(mesh->getFaces(i), quad, vol_discs); 
        surf_discs.push_back(surf_i);
      }
    }

    VolDiscPtr vol_disc;
    std::vector<SurfDiscPtr> surf_discs;
};

// creates mesh and a Discretization
class StandardDiscSetup : public StandardMeshBase
{
  public:
    virtual ~StandardDiscSetup() {}

    using StandardMeshBase::setup;

    virtual void setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& spec,
                       const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true}) override
    {
      StandardMeshBase::setup(quad_degree, sol_degree, spec, is_surf_dirichlet);
      disc = std::make_shared<Discretization>(mesh, quad_degree, quad_degree);
    }

  DiscPtr disc;
};

#endif
