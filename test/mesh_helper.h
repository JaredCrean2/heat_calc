#ifndef TEST_MESH_HELPER_H
#define TEST_MESH_HELPER_H

#include <memory>
//#include "mesh/mesh.h"
#include "mesh/mesh_input.h"
#include "discretization/surface_discretization.h"
#include "discretization/discretization.h"
#include "mesh/mesh_generator.h"
#include "quadrature.h"

namespace Mesh
{
  class MeshCG;
}

Mesh::MeshSpec getStandardMeshSpec();

std::vector<Mesh::MeshSpec> getStandardMeshSpecs();


std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const Mesh::MeshSpec& spec = getStandardMeshSpec(), int sol_degree=1,
                                const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true});

std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const std::vector<Mesh::MeshSpec>& meshspec, int sol_degree, 
                                               const std::vector<bool>& is_surf_dirichlet={true, true, true, true, true, true, true, true, true, true, true});


// creates a single block mesh with all dirichlet BCs
class StandardMeshBase
{
  public:
    virtual ~StandardMeshBase() {}

    // function must be called before first usage of any fields
    virtual void setup(const int quad_degree=5, int sol_degree=1, 
              const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true});


    virtual void setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& meshspec,
                       const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true});

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
                       const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true}) override;

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
                       const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true}) override;

    DiscPtr disc;
};

#include "mesh/mesh.h"

// creates mesh and a Discretization
class StandardDiscSetupMulti
{
  public:
    virtual ~StandardDiscSetupMulti() {}

    void setup(const int quad_degree, int sol_degree, const std::vector<Mesh::MeshSpec>& specs = getStandardMeshSpecs(),
               const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true, true, true, true, true, true});

    std::vector<Mesh::MeshSpec> specs;
    std::array<Real, 3> mesh_dim_mins;
    std::array<Real, 3> mesh_dim_maxs;
    std::shared_ptr<Mesh::MeshCG> mesh;
    Quadrature quad{getGaussianQuadrature(0)};

    std::vector<VolDiscPtr> vol_discs;
    std::vector<SurfDiscPtr> surf_discs;
    DiscPtr disc;
};

#endif
