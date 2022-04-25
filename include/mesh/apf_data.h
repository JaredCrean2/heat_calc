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

  explicit ApfData(apf::Mesh2* m, NumberingType* dof_nums=nullptr,
                   NumberingType* is_dirichlet=nullptr,
                   FSPtr sol_shape=nullptr,
                   FSPtr coord_shape=nullptr) :
    m(m), dof_nums(dof_nums), 
    is_dirichlet(is_dirichlet), 
    sol_shape(sol_shape.get()),
    coord_shape(coord_shape.get()),
    m_shared(makeSharedWithDeleter(m, apf::destroyMesh)),
    m_sol_shape(sol_shape),
    m_coord_shape(coord_shape)
  {
    apf::reorderMdsMesh(m);
  }

  ApfData(const ApfData&) = delete;

  ApfData& operator=(const ApfData&) = delete;

  apf::Mesh2* m;
  NumberingType* dof_nums       = nullptr;
  NumberingType* el_nums        = nullptr;
  NumberingType* is_dirichlet   = nullptr;  // if dof is dirichlet
  apf::FieldShape* sol_shape;    // FieldShape of solution
  apf::FieldShape* coord_shape;  // FieldShape of cooordinate field
  NumberingType* vol_groups     = nullptr;    // volume group numbering of elements

  // unfortunately, the apf API takes raw pointers rather than shared
  // pointers for everything, so we keep the raw pointer above, and also
  // have a shared pointer for memory management
  std::shared_ptr<apf::Mesh2> m_shared;
  FSPtr m_sol_shape;
  FSPtr m_coord_shape;
  REPtr m_ref_el_coord;
  REPtr m_ref_el_sol;
  std::shared_ptr<NumberingType> m_dof_nums;
  std::shared_ptr<NumberingType> m_el_nums;
  std::shared_ptr<NumberingType> m_is_dirichlet;
  std::shared_ptr<NumberingType> m_vol_groups;

  //std::vector<apf::MeshEntity*> elements;
};

}

#endif