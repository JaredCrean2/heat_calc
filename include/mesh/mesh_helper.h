#ifndef MESH_HELPER_H
#define MESH_HELPER_H

#include "mesh/mesh.h"

namespace Mesh {

void createGhostLayer(apf::Mesh2* m);


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

void setVolumeGroupNumbering(apf::Mesh2* m, const std::vector<MeshEntityGroupSpec>& vol_groups, ApfData::NumberingType* group_nums);

void getGroupElements(ApfData& apf_data, MeshEntityGroupSpec& volume_group,
                      std::vector<apf::MeshEntity*>& elements);

void getElementWeights(ApfData& apf_data, const std::vector<apf::MeshEntity*>& elements,
                        std::vector<int_least8_t>& element_weights);

void getGroupFaces(ApfData& apf_data, const MeshEntityGroupSpec& surf, REPtr ref_el_coord, 
                   const std::vector<MeshEntityGroupSpec>& volume_specs,
                   const std::vector<Index>& elnums_global_to_local,
                   std::vector<FaceSpec>& faces);

void getFaceWeights(ApfData& apf_data, const std::vector<FaceSpec>& faces, REPtr ref_el_coord,
                    std::vector<apf::MeshEntity*> elements,
                    std::vector<int_least8_t>& face_weights);

void getDofNums(ApfData& apf_data, const MeshEntityGroupSpec& vol_group,
                std::vector<apf::MeshEntity*>& elements_group,
                ArrayType<Index, 2>& nodenums);

void getDofNums(ApfData& apf_data, apf::MeshEntity* e, std::vector<int>& node_nums);


//const ArrayType<LocalIndex, 3>& getTensorProductMap(const int degree);

// get xi coordinates of coordinate nodes of 1D element
//const std::vector<Real>& getTensorProductXi(const int degree);

//const ArrayType<Real, 2>& getNormals();

void getCoords(ApfData& apf_data, const MeshEntityGroupSpec& vol_group,
               std::vector<apf::MeshEntity*>& elements_group,
               ArrayType<Real, 3>& coords);


void setDirichletDofs(ApfData& apf_data, int& dof_start);

// gives map relating nodes on the volume to the face .
// Output is nfaces (on reference element) x nnodes per face
//ArrayType<LocalIndex, 2> getFaceNodeMap(ApfData apf_data, apf::FieldShape* fshape);

//const ArrayType<LocalIndex, 2>& getFaceTensorProductMap(const int degree);

}
#endif
