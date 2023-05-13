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

void updateDependentDirichletValues(DiscVectorPtr disc_vec)
{
  auto dof_numbering = disc_vec->getDisc()->getDofNumbering();

  std::vector<Mesh::NodeTriplet> dest_nodes;
  for (int i=0; i < dof_numbering->getNumDirichletNodeSections(); ++i)
  {
    Mesh::NodeTriplet src_node = dof_numbering->getSrcDirichletNode(i);
    dof_numbering->getDestDirichletNodes(i, dest_nodes);

    Real src_val = disc_vec->getArray(src_node.vol_group)[src_node.el][src_node.node];
    for (auto& dest_node : dest_nodes)
    {
      disc_vec->getArray(dest_node.vol_group)[dest_node.el][dest_node.node] = src_val;
    }
  }

  // ordinarily disc_vec->markArrayModified() would be called here, but
  // the array entries modified do not appear in the vector
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