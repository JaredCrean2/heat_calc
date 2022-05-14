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
  int sol_degree              = 0;
  int coord_degree            = 0;
  int num_dofs                = 0;
  int num_dofs_total          = 0;
  int nodes_per_element       = 0;
  int coord_nodes_per_element = 0;
  int nodes_per_face          = 0;
};

class MeshCG
{
  using SInt = std::vector<VolumeGroup>::size_type;

  public:
    MeshCG(apf::Mesh2* m,
           std::vector<MeshEntityGroupSpec> volume_group_spec,
           std::vector<MeshEntityGroupSpec> bc_spec,
           std::vector<MeshEntityGroupSpec> other_surface_spec,
           const int solution_degree, const int coord_degree);

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
        if (m_all_face_spec[idx].getName() == name)
          return m_volume_spec[idx];

      throw std::invalid_argument(std::string("cannot find volume group named ") + name);
    }

    std::vector<apf::Vector3> normals_xi;

    //---------------------------------------------------------------------
    // dof functions

    // gets total number of dofs, including dirichlet BC dofs
    Index getNumTotalDofs() const {return m_dof_numbering.num_dofs_total;}

    // gets number of dofs, exluding dirichlet BC dofs
    Index getNumDofs() const {return m_dof_numbering.num_dofs;}

    // returns true if dof appears in the linear system (ie. not Dirichlet)
    bool isDofActive(const Index dof) const { return dof < getNumDofs();}

    void getElementDofs(const VolumeGroup& vol_group, int el_idx, std::vector<Index>& nodenums);

    // returns all the dofs connected to the given node (including self)
    void getDofConnectivity(const VolumeGroup& vol_group, const Index el_idx, const Index node, std::vector<DofInt>& dofs);

  private:
    void setVolumeIndices();

    void setSurfaceIndices(const std::vector<MeshEntityGroupSpec>& other_surface_spec);

    void setApfData();

    void createVolumeGroups();

    void createFaceGroups();


    // input
    ApfData m_apf_data;
    DofNumbering m_dof_numbering;
    std::vector<MeshEntityGroupSpec> m_volume_spec;
    std::vector<MeshEntityGroupSpec> m_bc_spec;
    std::vector<MeshEntityGroupSpec> m_all_face_spec;

    // computed data
    std::vector<VolumeGroup> m_vol_group;
    //std::vector<std::vector<FaceGroup>> m_bc_faces;
    std::vector<FaceGroup> m_all_faces;
    std::vector<apf::MeshEntity*> m_elements;
    std::vector<Index> m_elnums_global_to_local;  // element numbers to group element numbers
    std::vector< std::vector<Index> > m_elnums_local_to_global; // (group, local element number) -> global element number
    REPtr m_ref_el_coord;
    REPtr m_ref_el_sol;
    TensorProductMapper m_tensor_product_coord_map;
    TensorProductMapper m_tensor_product_sol_map;
};


} // namespace

#endif
