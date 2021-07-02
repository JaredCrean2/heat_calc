#ifndef TEST_MESH_HELPER_H
#define TEST_MESH_HELPER_H

#include <memory>
#include "mesh/mesh.h"
#include "mesh/mesh_input.h"

Mesh::MeshSpec getStandardMeshSpec();

std::shared_ptr<Mesh::MeshCG> makeStandardMesh(const Mesh::MeshSpec& spec = getStandardMeshSpec());

#endif
