#include "mesh_helper.h"
#include "mesh/mesh_generator.h"
#include "mesh/mesh_generator_multi.h"
#include "mesh/mesh.h"
#include "mesh/mesh_input.h"
#include "utils/error_handling.h"
#include <stdexcept>

//-----------------------------------------------------------------------------
// free functions

Mesh::MeshSpec getStandardMeshSpec()
{
  Mesh::MeshSpec meshspec = Mesh::getMeshSpec(0, 2, 0, 2, 0, 2, 5, 5, 5);

  return meshspec;
}

std::vector<Mesh::MeshSpec> getStandardMeshSpecs()
{
  Mesh::MeshSpec meshspec1 = Mesh::getMeshSpec(0, 2, 0, 2, 0, 2, 2, 2, 2);
  Mesh::MeshSpec meshspec2 = Mesh::getMeshSpec(0, 2, 2, 4, 0, 2, 2, 3, 2);

  return {meshspec1, meshspec2};
}


std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const Mesh::MeshSpec& meshspec, int sol_degree, const std::vector<bool>& is_surf_dirichlet)
{
  assertAlways(is_surf_dirichlet.size() == 6, "is_surf_dirichlet must have length 6");

  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  auto m = Mesh::make_parallel_mesh(meshspec, comm_size, &(Mesh::identity));
  //auto generator = Mesh::make_mesh_generator(meshspec, &(Mesh::identity));
  //auto m = generator.generate();

  // make volume groups
  Mesh::MeshEntityGroupSpec volume_group("volume0");
  volume_group.addModelEntity(Mesh::ModelEntitySpec(3, 0));
  std::vector<Mesh::MeshEntityGroupSpec> volume_groups{volume_group};

  // make surface groups
  std::vector<Mesh::MeshEntityGroupSpec> surface_groups;
  for (int i=0; i < 6; ++i)
  {
    surface_groups.emplace_back(std::string("surface") + std::to_string(i));
    surface_groups.back().addModelEntity(Mesh::ModelEntitySpec(2, i), Mesh::ModelEntitySpec(3, 0));
    surface_groups.back().setIsDirichlet(is_surf_dirichlet[i]);
  }

  std::vector<Mesh::MeshEntityGroupSpec> other_surfaces;
  auto mesh = Mesh::createMeshCG(m, volume_groups, surface_groups,
                                 other_surfaces, sol_degree, 1);

  return mesh;
}

std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const std::vector<Mesh::MeshSpec>& meshspec, int sol_degree, const std::vector<bool>& is_surf_dirichlet)
{
  assertAlways(is_surf_dirichlet.size() == 11, "is_surf_dirichlet must have length 11");

  //auto generator = Mesh::make_mesh_generator(meshspec, &(Mesh::identity));
  //auto m = generator.generate();
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  auto m = Mesh::make_parallel_mesh(meshspec, comm_size, &(Mesh::identity));

  // make volume groups
  std::vector<Mesh::MeshEntityGroupSpec> volume_groups;
  for (int i=0; i < 2; ++i) 
  {
    Mesh::MeshEntityGroupSpec volume_group(std::string("volume") + std::to_string(i));
    volume_group.addModelEntity(Mesh::ModelEntitySpec(3, i));
    volume_groups.push_back(volume_group);
  }

  // make surface groups
  std::vector<Mesh::MeshEntityGroupSpec> surface_groups, internal_groups;
  int idx = 0;
  for (int i=0; i < 11; ++i)
  {
    if (i == 2)
      internal_groups.emplace_back(std::string("surface") + std::to_string(i));
    else
    {
      surface_groups.emplace_back(std::string("surface") + std::to_string(i));
      surface_groups.back().setIsDirichlet(is_surf_dirichlet[idx]);
      idx++;
    }
    int vol_group = i < 6 ? 0 : 1;

    auto& group = i == 2 ? internal_groups.back() : surface_groups.back();
    group.addModelEntity(Mesh::ModelEntitySpec(2, i), Mesh::ModelEntitySpec(3, vol_group));
  }

  auto mesh = Mesh::createMeshCG(m, volume_groups, surface_groups,
                                 internal_groups, sol_degree, 1);

  return mesh;
}


// StandardMeshBase
//-----------------------------------------------------------------------------
void StandardMeshBase::setup(const int quad_degree, int sol_degree, 
           const std::vector<bool>& is_surf_dirichlet)
{
  spec = getStandardMeshSpec();
  spec.nx = 4; spec.ny = 5, spec.nz = 6;

  //spec.nx = 2; spec.ny = 1, spec.nz = 1;
  //spec.xmin = 0; spec.ymin = 0; spec.zmin = 0;
  //spec.xmax = 1; spec.ymax = 1; spec.zmax = 1;
  setup(quad_degree, sol_degree, spec, is_surf_dirichlet);
}


void StandardMeshBase::setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& meshspec,
          const std::vector<bool>& is_surf_dirichlet)
{
  spec = meshspec;
  //quad = getGaussianQuadrature(quad_degree);

  mesh_dim_mins = {spec.xmin, spec.ymin, spec.zmin};
  mesh_dim_maxs = {spec.xmax, spec.ymax, spec.zmax};
  mesh = makeStandardMesh(spec, sol_degree, is_surf_dirichlet);
}


// StandardMeshSetup
//-----------------------------------------------------------------------------
void StandardMeshSetup::setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& spec,
                              const std::vector<bool>& is_surf_dirichlet)
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

// StandardDiscSetup
//-----------------------------------------------------------------------------
void StandardDiscSetup::setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& spec,
                              const std::vector<bool>& is_surf_dirichlet)
{
  StandardMeshBase::setup(quad_degree, sol_degree, spec, is_surf_dirichlet);
  disc = std::make_shared<Discretization>(mesh, quad_degree, quad_degree);
}

// StandardDiscSetupMulti
//-----------------------------------------------------------------------------

void StandardDiscSetupMulti::setup(const int quad_degree, int sol_degree, const std::vector<Mesh::MeshSpec>& specs,
                                   const std::vector<bool>& is_surf_dirichlet)
{
  this->specs = specs;
  quad = getGaussianQuadrature(quad_degree);

  mesh_dim_mins = {specs[0].xmin, specs[0].ymin, specs[0].zmin};
  mesh_dim_maxs = {specs[0].xmax, specs.back().ymax, specs[0].zmax};
  mesh = makeStandardMesh(specs, sol_degree, is_surf_dirichlet);

  disc = std::make_shared<Discretization>(mesh, quad_degree, quad_degree);

  vol_discs.resize(0);
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
    vol_discs.push_back(disc->getVolDisc(i));

  surf_discs.resize(0); 
  for (int i=0; i < disc->getNumSurfDiscs(); ++i)
    surf_discs.push_back(disc->getSurfDisc(i));
}