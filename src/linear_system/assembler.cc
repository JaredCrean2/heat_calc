#include "linear_system/assembler.h"
#include "mesh/mesh.h"

namespace linear_system {

static_assert(std::is_signed<DofInt>::value, "DofInt must be a signed integral type");

Assembler::Assembler(DiscPtr disc, LargeMatrixPtr mat) :
  m_disc(disc),
  m_dof_nums(disc->getDofNumbering()),
  m_alpha(1),
  m_vol_dofs(disc->getVolDisc(0)->getNumSolPtsPerElement()),
  m_face_dofs(disc->getSurfDisc(0)->getNumSolPtsPerFace()),
  m_matrix(mat)
{
  disc->getMesh()->getLocalToGlobalDofs(m_local_dof_to_global);
  
}


void Assembler::assembleVolume(int vol_disc_idx, int elnum, ArrayType<Real, 2>& jac)
{

  auto vol_disc = m_disc->getVolDisc(vol_disc_idx);
  assert(jac.shape()[0] == vol_disc->getNumSolPtsPerElement());
  assert(jac.shape()[1] == vol_disc->getNumSolPtsPerElement());

  // TODO: this could be a contiguous view of dofs.  The only reason it isn't
  //       is assembleVolume is a virtual method, and methods can't be both
  //       template and virtual
  auto& dofs = m_dof_nums->getDofs(vol_disc_idx);
  for (unsigned int i=0; i < dofs.shape()[1]; ++i)
  {
    auto dof = dofs[elnum][i];
    if (!m_dof_nums->isDofActive(dof))
      dof = -1;
    else
      dof = m_local_dof_to_global[dof];

    m_vol_dofs[i] = dof;

    for (unsigned int j=0; j < dofs.shape()[1]; ++j)
      jac[i][j] *= m_alpha;
  }

  m_matrix->assembleValues(m_vol_dofs, jac);
}

void Assembler::assembleFace(int surf_disc_idx, int facenum, ArrayType<Real, 2>& jac)
{
  auto surf_disc = m_disc->getSurfDisc(surf_disc_idx);
  assert(jac.shape()[0] == surf_disc->getNumSolPtsPerFace());
  assert(jac.shape()[1] == surf_disc->getNumSolPtsPerFace());

  auto& face_nodemap = surf_disc->face_group.getFaceNodesSol();
  auto& face_spec    = surf_disc->face_group.faces[facenum];
  auto& dofs         = m_dof_nums->getDofs(face_spec.vol_group);
  for (int i=0; i < surf_disc->getNumSolPtsPerFace(); ++i)
  {
    auto dof = dofs[face_spec.el_group][ face_nodemap[face_spec.face][i] ];
    if (m_dof_nums->isDofActive(dof))
      m_face_dofs[i] = m_local_dof_to_global[dof];
    else
      m_face_dofs[i] = -1;

    for (int j=0; j < surf_disc->getNumSolPtsPerFace(); ++j)
      jac[i][j] *= m_alpha;
  }

  m_matrix->assembleValues(m_face_dofs, jac);
}

}  // namespace