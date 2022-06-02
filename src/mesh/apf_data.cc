#include "mesh/apf_data.h"
#include "mesh/reference_element_apf.h"

namespace Mesh {

#ifdef MESH_USE_MDS_NUMBERING
apf::ApfMDSNumbering* createNumbering(apf::Mesh2* mesh, const char* name, apf::FieldShape* shape, int components)
{
  return apf::createNumberingMDS(mesh, name, shape, components);
}
#else
apf::Numbering* createNumbering(apf::Mesh2* mesh, const char* name, apf::FieldShape* shape, int components)
{
  return apf::createNumbering(mesh, name, shape, components);
}
#endif



ApfData::ApfData(apf::Mesh2* m,
                REPtr ref_el_sol,
                REPtr ref_el_coord) :

  m_shared(makeSharedWithDeleter(m, apf::destroyMesh)),
  m_sol_shape(apf::getHexFieldShape(ref_el_sol)),
  m_coord_shape(apf::getHexFieldShape(ref_el_coord)),
  m_ref_el_coord(ref_el_coord),
  m_ref_el_sol(ref_el_sol),

  m(m),
  sol_shape(m_sol_shape.get()),
  coord_shape(m_coord_shape.get()),

  dof_nums(createNumbering(m, "local_dof_nums",  sol_shape, 1)),
  global_dof_nums(createNumbering(m, "global_dof_nums", sol_shape, 1)),
  el_nums(createNumbering(m, "el_nums", apf::getConstant(3), 1)),
  is_dirichlet(createNumbering(m, "is_dirichlet", this->sol_shape, 1)),
  vol_groups(createNumbering(m, "vol_group", apf::getConstant(3), 1))
{
  using DeleteNumbering = void (*)(ApfData::NumberingType*);
  DeleteNumbering deleter      = &apf::destroyNumbering;
  m_dof_nums        = makeSharedWithDeleter(dof_nums, deleter);
  m_global_dof_nums = makeSharedWithDeleter(global_dof_nums, deleter);
  m_el_nums         = makeSharedWithDeleter(el_nums, deleter);
  m_is_dirichlet    = makeSharedWithDeleter(is_dirichlet, deleter);
  m_vol_groups      = makeSharedWithDeleter(vol_groups, deleter);
}

}