#ifndef MESH_HELPER_H
#define MESH_HELPER_H

#include "mesh/mesh.h"

namespace Mesh {

// Handling Dirichlet BCs

// Populates apf_data.is_dirichlet with true for all dofs on Dirichlet faces
void setDirichletNodes(ApfData& apf_data, std::vector<MeshEntityGroupSpec>& m_bc_spec);

  // get the ModelEntitySpec for a given MeshEntity
ModelEntitySpec getMESpec(apf::Mesh* m, apf::MeshEntity* e);

// sets apf_data.is_dirichlet  all dofs affecting a given MeshEntity to
// the given value
void setDofsDirichlet(ApfData& apf_data, apf::MeshEntity* e, const int val);


//TODO: need to do both local and global dof numbering
//TODO: need to setup parallel communication

//-----------------------------------------------------------------------------
// Dof numbering

//void initializeElnums(ApfData& apf_data);

// counts number of elements in a volume group
// Note: this iterates over all elements in the mesh, and may be slow
int countNumEls(ApfData& apf_data, const MeshEntityGroupSpec& vol_group);

void setVolumeGroupNumbering(apf::Mesh2* m, const std::vector<MeshEntityGroupSpec>& vol_groups, apf::Numbering* group_nums);

//TODO: implement this
// on exit, elnums_start and elnums_end are the starting numbers for the next
// block
void setDofNums(apf::Mesh* m, MeshEntityGroupSpec volume_group,
                   int& elnums_start, int& dof_start,
                   apf::Numbering* is_dirichlet, apf::Numbering* dof_nums,
                   apf::Numbering* el_nums);

void getDofNums(ApfData apf_data, const MeshEntityGroupSpec& vol_group,
                ArrayType<Index, 2>& nodenums,
                std::vector<apf::MeshEntity*>& elements_group);


void setDirichletDofs(ApfData& apf_data, int& dof_start);


ArrayType<LocalIndex, 2> getFaceNodeMap(ApfData apf_data);

}
#endif
