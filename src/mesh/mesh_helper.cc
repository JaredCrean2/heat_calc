#include "mesh/mesh_helper.h"
#include "apfShape.h"
#include "apf.h"
#include <iostream>

namespace Mesh {

//-----------------------------------------------------------------------------
// Handling Dirichlet BCs

// Populates apf_data.is_dirichlet with true for all dofs on Dirichlet faces
void setDirichletNodes(ApfData& apf_data, std::vector<MeshEntityGroupSpec>& m_bc_spec)
{
  // initialize is_dirichlet to false
  for (int dim=0; dim <= 3; ++dim)
    if (apf_data.sol_shape->hasNodesIn(dim))
    {
      apf::MeshIterator* it = apf_data.m->begin(dim);
      while ( apf::MeshEntity* e = apf_data.m->iterate(it) )
      {
        auto type = apf_data.m->getType(e);
        for (int i=0; i < apf_data.sol_shape->countNodesOn(type); ++i)
          apf::number(apf_data.is_dirichlet, e, i, 0, false);

      }  
      apf_data.m->end(it);
    } 

  // set is_dirichlet to true for dirichlet surfaces
  for (const auto& surf : m_bc_spec)
  {
    if (!surf.getIsDirichlet())
      continue;

    apf::MeshIterator* it = apf_data.m->begin(2);
    apf::MeshEntity* e;
    while ( (e = apf_data.m->iterate(it)) )
    {
      auto me_spec = getMESpec(apf_data.m, e);
      if (surf.hasModelEntity(me_spec))
        setDofsDirichlet(apf_data, e, true);
    }
    apf_data.m->end(it);
  }
}

// get the ModelEntitySpec for a given MeshEntity
ModelEntitySpec getMESpec(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::ModelEntity* me = m->toModel(e);
  return ModelEntitySpec(m->getModelType(me), m->getModelTag(me));
}


// sets apf_data.is_dirichlet  all dofs affecting a given MeshEntity to
// the given value
void setDofsDirichlet(ApfData& apf_data, apf::MeshEntity* e, const int val)
{
  apf::Downward down;
  int ndown;
  int e_dim = apf::Mesh::typeDimension[apf_data.m->getType(e)];
  for (int dim = 0; dim <= e_dim; ++dim)
    if (apf_data.sol_shape->hasNodesIn(dim))
    {
      // TODO: does Pumi define entities to be downward adjacent to themselves?
      if (dim == e_dim)
      {
        ndown = 1;
        down[0] = e;
      } else
        ndown = apf_data.m->getDownward(e, dim, down);

      for (int i=0; i < ndown; ++i)
      {
        auto type = apf_data.m->getType(down[i]);
        for (int j=0; j < apf_data.sol_shape->countNodesOn(type); ++j)
          apf::number(apf_data.is_dirichlet, down[i], j, 0, true);
      }
    }
}

//TODO: need to do both local and global dof numbering
//TODO: need to setup parallel communication

//-----------------------------------------------------------------------------
// Dof numbering
/*
void initializeElnums(ApfData& apf_data)
{

  std::size_t val = apf_data.m->count(3);
  apf::MeshIterator* it = apf_data.m->begin(3);
  apf::MeshEntity* e;
  while( (e = apf_data.m->iterate(it)) )
    apf::number(apf_data.el_nums, e, 0, 0, val);

  apf_data.m->end(it);
}
*/

int countNumEls(ApfData& apf_data, const MeshEntityGroupSpec& vol_group)
{
  int numel = 0;
  apf::MeshIterator* it = apf_data.m->begin(3);
  apf::MeshEntity* e;
  while ( (e = apf_data.m->iterate(it)) )
  {
    auto me = getMESpec(apf_data.m, e);
    if (vol_group.hasModelEntity(me))
      numel += 1;
  }

  return numel;
}

void setVolumeGroupNumbering(apf::Mesh2* m, const std::vector<MeshEntityGroupSpec>& vol_groups, apf::Numbering* group_nums)
{
  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;

  while ( (e = m->iterate(it)) )
  {
    auto me = getMESpec(m, e);
    for (const auto& vol_group : vol_groups)
      if (vol_group.hasModelEntity(me))
      {
        apf::number(group_nums, e, 0, 0, vol_group.getIdx());
        break;
      }
  }
}

void getGroupElements(ApfData apf_data, MeshEntityGroupSpec& volume_group,
                      std::vector<apf::MeshEntity*>& elements)
{
  apf::MeshIterator* it = apf_data.m->begin(3);
  apf::MeshEntity* e;
  while ( (e = apf_data.m->iterate(it)) )
  {
    auto me_spec = getMESpec(apf_data.m, e);
    if (!(volume_group.hasModelEntity(me_spec)))
      continue;

    elements.push_back(e);
  }
  apf_data.m->end(it);
}

void getDofNums(ApfData apf_data, const MeshEntityGroupSpec& vol_group,
                std::vector<apf::MeshEntity*>& elements,
                ArrayType<Index, 2>& nodenums)
{

  apf::Downward down;
  int el_idx = 0;
  for( auto e : elements )
  {
    int idx=0, offset=0;
    for (int dim=0; dim <= 3; ++dim)
    {
      int ndown = apf_data.m->getDownward(e, dim, down);

      int nnodes_dim = apf_data.sol_shape->countNodesOn(apf_data.m->getType(down[0]));
      for (int i=0; i < ndown; ++i)
        for (int j=0; j < nnodes_dim; ++j)
        {
          nodenums[el_idx][idx++] =
            apf::getNumber(apf_data.dof_nums, down[i], j, 0);
        }

      offset += ndown * nnodes_dim;
    }
    el_idx += 1;
  }
}


void getCoords(ApfData apf_data, const MeshEntityGroupSpec& vol_group,
               std::vector<apf::MeshEntity*>& elements,
               ArrayType<Real, 3>& coords)
{
  apf::Downward down;
  apf::FieldShape* fshape = apf_data.m->getShape();
  apf::Vector3 point;
  int el_idx = 0;
  for( auto e : elements )
  {
    int idx=0, offset=0;
    for (int dim=0; dim <= 3; ++dim)
    {
      int ndown = apf_data.m->getDownward(e, dim, down);

      int nnodes_dim = fshape->countNodesOn(apf_data.m->getType(down[0]));
      for (int i=0; i < ndown; ++i)
        for (int j=0; j < nnodes_dim; ++j)
        {
          apf_data.m->getPoint(down[i], j, point);

          coords[el_idx][idx][0] = point.x();
          coords[el_idx][idx][1] = point.y();
          coords[el_idx][idx][2] = point.z();
          ++idx;
        }

      offset += ndown * nnodes_dim;
    }
    ++el_idx;
  }
}



void setDirichletDofs(ApfData& apf_data, int& dof_start)
{
  for (int dim=0; dim <=3; ++dim)
    if (apf_data.sol_shape->hasNodesIn(dim))
    {
      apf::MeshIterator* it = apf_data.m->begin(dim);
      apf::MeshEntity* e;
      while ( (e = apf_data.m->iterate(it)) )
      {
        int type = apf_data.m->getType(e);
        for (int i=0; i < apf_data.sol_shape->countNodesOn(type); ++i)
          if (apf::getNumber(apf_data.is_dirichlet, e, i, 0))
            apf::number(apf_data.dof_nums, e, i, 0, dof_start++);
      }
    }
}


ArrayType<LocalIndex, 2> getFaceNodeMap(ApfData apf_data)
{
  apf::Downward faces;
  apf::Downward face_entities;
  apf::Downward el_entities;

  apf::MeshIterator* it = apf_data.m->begin(3);
  apf::MeshEntity* e = apf_data.m->iterate(it);
  apf_data.m->end(it);

  int nfaces = apf_data.m->getDownward(e, 2, faces);
  auto type = apf_data.m->getType(faces[0]);
  apf::EntityShape* eshape = apf_data.sol_shape->getEntityShape(type);
  int nnodes = eshape->countNodes();


  // note: this must be consistent with how getDofs maps dofs to a 1D array
  // The convention here is: vertices, then edges, then faces (in APF
  // topology order)

  ArrayType<LocalIndex, 2> nodemap(boost::extents[nfaces][nnodes]);

  for (int face=0; face < nfaces; ++face)
  {
    int idx = 0, offset = 0;
    for (int dim=0; dim <= 2; ++dim)
    {
      int n_el_entities, n_face_entities;
      if (dim == 2)
      {
        n_el_entities = nfaces;
        for (int i=0; i < nfaces; ++i)
          el_entities[i] = faces[i];

        n_face_entities = 1;
        face_entities[0] = faces[face];
      } else
      { 
        n_el_entities = apf_data.m->getDownward(e, dim, el_entities);
        n_face_entities = apf_data.m->getDownward(faces[face], dim, face_entities);
      }

      int nnodes_dim = apf_data.sol_shape->countNodesOn(apf_data.m->getType(face_entities[0]));
      for (int i=0; i < n_face_entities; ++i)
      {
        int foundidx = 0;
        for (int j=0; j < n_el_entities; ++j)
          if (face_entities[i] == el_entities[j])
          {
            foundidx = j;
            break;
          }

        for (int node=0; node < nnodes_dim; ++node)
          nodemap[face][idx++] = foundidx * nnodes_dim + node + offset;
      }
      offset += n_el_entities * nnodes_dim;
    }
  }

  return nodemap;
}

} // namespace
