#ifndef MESH_H
#define MESH_H

#include "ProjectDefs.h"
#include "utils/memory.h"
#include "mesh/mesh_input.h"
#include "mesh/reference_element.h"  //TODO: remove this
#include "mesh/reference_element_interface.h"
#include "mesh/reference_element_apf.h"
#include "mesh/tensor_product_mapper.h"
#include "mesh/volume_group.h"
#include "mesh/face_group.h"
#include "mesh/apf_data.h"
#include "mesh/apfMDSField.h"
#include "mesh/FieldDataManager.h"

#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfMDS.h>
#include <vector>
#include <iostream>


namespace Mesh {

using REPtr = reference_element::REPtr;

//TODO: merge into main initialize method
bool initialize();

struct DofNumbering
{
  int sol_degree                        = 0;
  int coord_degree                      = 0;
  int num_dofs_local                    = 0;
  int num_dofs_owned                    = 0;
  int num_dofs_total                    = 0;  // local + dirichlet dofs
  int local_owned_to_global_dof_offset  = 0;
  int nodes_per_element                 = 0;
  int coord_nodes_per_element           = 0;
  int nodes_per_face                    = 0;
};

class MeshCG
{
  using SInt = std::vector<VolumeGroup>::size_type;

    MeshCG(apf::Mesh2* m,
           std::vector<MeshEntityGroupSpec> volume_group_spec,
           std::vector<MeshEntityGroupSpec> bc_spec,
           std::vector<MeshEntityGroupSpec> other_surface_spec,
           const int solution_degree, const int coord_degree);

  public:
    // getting faces
    const FaceGroup& getFaces(const MeshEntityGroupSpec& surf)
    {
      return m_all_faces.at(surf.getIdx());
    }

    const FaceGroup& getFaces(const Index idx)
    {
      return m_all_faces.at(idx);
    }

    Index getNumBCSurfaces() const {return m_bc_spec.size();}

    Index getNumSurfaces() const {return m_all_face_spec.size();}

    const MeshEntityGroupSpec& getSurface(const std::string& name) const
    {
      for (SInt idx=0; idx < m_all_face_spec.size(); ++idx)
        if (m_all_face_spec[idx].getName() == name)
          return m_all_face_spec[idx];

      throw std::invalid_argument(std::string("cannot find surface named ") + name);
    }

    // getting elements
    VolumeGroup& getElements(const MeshEntityGroupSpec& surf)
    {
      return m_vol_group.at(surf.getIdx());
    }

    VolumeGroup& getElements(const Index idx)
    {
      return m_vol_group.at(idx);
    }

    Index getNumVolumeGroups() const {return m_vol_group.size();}

    const MeshEntityGroupSpec& getVolumeGroup(const std::string& name) const
    {
      for (SInt idx=0; idx < m_volume_spec.size(); ++idx)
      {
        if (m_vol_group[idx].getName() == name)
          return m_volume_spec[idx];
      }

      throw std::invalid_argument(std::string("cannot find volume group named ") + name);
    }

    std::vector<apf::Vector3> normals_xi;

    //---------------------------------------------------------------------
    // dof functions

    // gets total number of local dofs, including dirichlet BC dofs
    Index getNumTotalDofs() const {return m_dof_numbering.num_dofs_total;}

    // gets number of local dofs, exluding dirichlet BC dofs
    Index getNumDofs() const {return m_dof_numbering.num_dofs_local;}

    // get number of owned dofs, excluding dirichlet BC dofs
    Index getNumOwnedDofs() const { return m_dof_numbering.num_dofs_owned; }

    // returns true if dof appears in the linear system (ie. not Dirichlet)
    // dof must be a *local* dof
    bool isDofActive(const Index dof) const
    { 
      assert(dof >= 0 && dof < getNumTotalDofs());
      return dof < getNumDofs();
    }

    // returns the local dofs on a given element
    void getElementDofs(const VolumeGroup& vol_group, int el_idx, std::vector<Index>& nodenums);

    // returns all the local dofs connected to the given node (including self)
    void getDofConnectivity(const VolumeGroup& vol_group, const Index el_idx, const Index node, std::vector<DofInt>& dofs);

    // for non-owned dofs, gives both the global dof number and the corresponding local dof number
    void getGhostDofInfo(std::vector<DofInt>& global_dofs, std::vector<DofInt>& local_dofs);
    
    // returns a vector v such that v[owned_dof_num] = local_dof_num
    void getOwnedLocalDofInfo(std::vector<DofInt>& owned_local_dofs);

    void getLocalToGlobalDofs(std::vector<DofInt>& local_to_global_dofs);

    FieldDataManager& getFieldDataManager() { return m_field_data_manager; }

    void writeVtkFiles(const std::string& fname);

  private:
    void setVolumeIndices();

    void setSurfaceIndices(const std::vector<MeshEntityGroupSpec>& other_surface_spec);

    void setApfData();

    void createVolumeGroups();

    void createFaceGroups();


    // input
    std::vector<MeshEntityGroupSpec> m_volume_spec;
    std::vector<MeshEntityGroupSpec> m_bc_spec;
    std::vector<MeshEntityGroupSpec> m_all_face_spec;
    REPtr m_ref_el_coord;
    REPtr m_ref_el_sol;
    ApfData m_apf_data;


    // computed data
    DofNumbering m_dof_numbering;
    std::vector<VolumeGroup> m_vol_group;
    //std::vector<std::vector<FaceGroup>> m_bc_faces;
    std::vector<FaceGroup> m_all_faces;
    std::vector<apf::MeshEntity*> m_elements;
    std::vector<Index> m_elnums_global_to_local;  // element numbers to group element numbers
    std::vector< std::vector<Index> > m_elnums_local_to_global; // (group, local element number) -> global element number

    TensorProductMapper m_tensor_product_coord_map;
    TensorProductMapper m_tensor_product_sol_map;
    FieldDataManager m_field_data_manager;

    friend std::shared_ptr<MeshCG> createMeshCG(apf::Mesh2* m,
                                     std::vector<MeshEntityGroupSpec> volume_group_spec,
                                     std::vector<MeshEntityGroupSpec> bc_spec,
                                     std::vector<MeshEntityGroupSpec> other_surface_spec,
                                     const int solution_degree, const int coord_degree);
};

std::shared_ptr<MeshCG> createMeshCG(apf::Mesh2* m,
                                     std::vector<MeshEntityGroupSpec> volume_group_spec,
                                     std::vector<MeshEntityGroupSpec> bc_spec,
                                     std::vector<MeshEntityGroupSpec> other_surface_spec,
                                     const int solution_degree, const int coord_degree);

} // namespace

#endif
