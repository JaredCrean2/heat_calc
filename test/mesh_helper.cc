#include "mesh_helper.h"
#include "mesh/mesh_generator.h"
#include "mesh/mesh.h"
#include "mesh/mesh_input.h"
#include "utils/error_handling.h"
#include <stdexcept>

Mesh::MeshSpec getStandardMeshSpec()
{
  Mesh::MeshSpec meshspec = Mesh::getMeshSpec(0, 2, 0, 2, 0, 2, 5, 5, 5);

  return meshspec;
}


std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const Mesh::MeshSpec& meshspec, int sol_degree, const std::vector<bool>& is_surf_dirichlet)
{
  assertAlways(is_surf_dirichlet.size() == 6, "is_surf_dirichlet must have length 6");

  auto generator = Mesh::make_mesh_generator(meshspec, &(Mesh::identity));
  auto m = generator.generate();

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
  auto mesh = std::make_shared<Mesh::MeshCG>(m, volume_groups, surface_groups,
                                       other_surfaces, sol_degree, 1);

  return mesh;
}
