#ifndef MESH_APF_DATA_H
#define MESH_APF_DATA_H

#include "apf.h"
#include "apfNumbering.h"
#include "mesh/reference_element_interface.h"
#include "mesh/reference_element_apf.h"
#include "mesh/apfMDSField.h"
#include "utils/memory.h"


#include <memory>


namespace Mesh {

using REPtr = reference_element::REPtr;




struct ApfData
{
  using FSPtr = std::shared_ptr<apf::FieldShapeRefEl>;
#ifdef MESH_USE_MDS_NUMBERING
  using NumberingType = apf::ApfMDSNumbering;
#else
  using NumberingType = apf::Numbering;
#endif

  explicit ApfData(apf::Mesh2* m,
                   REPtr ref_el_sol,
                   REPtr ref_el_coord);

  ApfData(const ApfData&) = delete;

  ApfData& operator=(const ApfData&) = delete;

  // unfortunately, the apf API takes raw pointers rather than shared
  // pointers for everything, so we keep the raw pointer and also
  // have a shared pointer for memory management
  std::shared_ptr<apf::Mesh2> m_shared;
  FSPtr m_sol_shape;
  FSPtr m_coord_shape;
  REPtr m_ref_el_coord;
  REPtr m_ref_el_sol;

  apf::Mesh2* m;
  apf::FieldShape* sol_shape;    // FieldShape of solution
  apf::FieldShape* coord_shape;  // FieldShape of cooordinate field

  NumberingType* dof_nums        = nullptr;
  NumberingType* global_dof_nums = nullptr;
  NumberingType* el_nums         = nullptr;
  NumberingType* is_dirichlet    = nullptr;  // if dof is dirichlet
  NumberingType* vol_groups      = nullptr;    // volume group numbering of elements


  std::shared_ptr<NumberingType> m_dof_nums;
  std::shared_ptr<NumberingType> m_global_dof_nums;
  std::shared_ptr<NumberingType> m_el_nums;
  std::shared_ptr<NumberingType> m_is_dirichlet;
  std::shared_ptr<NumberingType> m_vol_groups;





  //std::vector<apf::MeshEntity*> elements;
};

}

#endif