#include "discretization/DirichletBC.h"
#include "mesh/mesh.h"


void applyDirichletValues(DirichletBCPtr bc, const Real t, DiscVectorPtr disc_vec)
{
  //std::cout << "\nEntered applyDirichletValues" << std::endl;
  auto surf = bc->getSurfDisc();
  std::vector<Real> vals(surf->getNumSolPtsPerFace());
  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    bc->getValue(face, t, vals.data());
    auto& face_spec = surf->face_group.faces[face];

    auto& arr = disc_vec->getArray(face_spec.vol_group);
    auto& nodemap = surf->face_group.getFaceNodesSol();
    for (int j=0; j < surf->getNumSolPtsPerFace(); ++j)
    {
      arr[face_spec.el_group][nodemap[face_spec.face][j]] = vals[j];
    }
  }
}

namespace impl {

Real errorFunc(Real x, Real y, Real z, Real t)
{
  throw std::runtime_error("unsteady Dirichlet BC not set up");
}

Real zeroFunc(Real x, Real y, Real z, Real t)
{
  return 0;
}

}