#include "mesh/mesh_helper.h"
#include "mesh/apfMDSField.h"
#include "mesh/apfMDSFieldReduction.h"
#include "apfShape.h"
#include "apf.h"
#include "pumi.h"
#include "mesh/mesh_input.h"
#include <iostream>

namespace Mesh {

void createGhostLayer(apf::Mesh2* m)
{
  std::shared_ptr<gModel> pumi_geom_shared = std::make_shared<gModel>(m->getModel());

  pumi::instance()->mesh = m;
  pumi::instance()->model = pumi_geom_shared.get();
  pumi_mesh_setCount(m, nullptr);
  int num_layers = 1;
  pumi_ghost_createLayer(m, 0, 3, num_layers, true);
  apf::reorderMdsMesh(m);
}


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

  auto reducer = fast_field::make_field_reduction(*(apf_data.is_dirichlet), std::logical_or<int>());
  reducer.execute();
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
          apf::number(apf_data.is_dirichlet, down[i], j, 0, val);
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

  apf_data.m->end(it);

  return numel;
}

void setVolumeGroupNumbering(apf::Mesh2* m, const std::vector<MeshEntityGroupSpec>& vol_groups, ApfData::NumberingType* group_nums)
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
  m->end(it);
}

void getGroupElements(ApfData& apf_data, MeshEntityGroupSpec& volume_group,
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

void getElementWeights(ApfData& apf_data, const std::vector<apf::MeshEntity*>& elements,
                        std::vector<int_least8_t>& element_weights)
{
  element_weights.resize(elements.size());
  for (size_t i=0; i < elements.size(); ++i)
    element_weights[i] = apf_data.m->isGhost(elements[i]) ? 0 : 1;
}


void getGroupFaces(ApfData& apf_data, const MeshEntityGroupSpec& surf, REPtr ref_el_coord, 
                   const std::vector<MeshEntityGroupSpec>& volume_specs,
                   const std::vector<Index>& elnums_global_to_local,
                   std::vector<FaceSpec>& faces)
{
  apf::Downward down;
  apf::MeshIterator* it = apf_data.m->begin(2);
  apf::MeshEntity* e;
  while ( (e = apf_data.m->iterate(it)) )
  {
    auto me = getMESpec(apf_data.m, e);
    if (surf.hasModelEntity(me))
    {
      auto me_parent = surf.getParentEntity(me);
      for (int i=0; i < apf_data.m->countUpward(e); ++i)
      {
        apf::MeshEntity* e_i = apf_data.m->getUpward(e, i);
        auto me_i = getMESpec(apf_data.m, e_i);
        if (me_i == me_parent)
        {
          // get local face number
          int nfaces_i = apf_data.m->getDownward(e_i, 2, down);
          int localidx = -1;
          for (int j=0; j < nfaces_i; ++j)
            if (down[j] == e)
            {
              localidx = j;
              break;
            }

          assert(localidx != -1);
          localidx = ref_el_coord->getREEntityIndex(2, localidx);

          int elnum = apf::getNumber(apf_data.el_nums, e_i, 0, 0);
          int elnum_local = elnums_global_to_local[elnum];
          int vol_group_idx = -1;
          // TODO: could cache this, because it is a function of me_parent
          for (const auto& vol_group : volume_specs)
            if (vol_group.hasModelEntity(me_parent))
              vol_group_idx = vol_group.getIdx();

          assert(vol_group_idx != -1);

          faces.push_back(FaceSpec(elnum, elnum_local, localidx, vol_group_idx));
          break;
        }
      }
    }
  }
  apf_data.m->end(it);
}

void getFaceWeights(ApfData& apf_data, const std::vector<FaceSpec>& faces, REPtr ref_el_coord,
                    std::vector<apf::MeshEntity*> elements,
                    std::vector<int_least8_t>& face_weights)
{
  apf::Downward down;
  face_weights.resize(faces.size());
  for (size_t i=0; i < faces.size(); ++i)
  {
    int local_face_idx = faces[i].face;
    int local_face_apf = ref_el_coord->getApfEntityIndex(2, local_face_idx);
    apf::MeshEntity* el_apf = elements[faces[i].el];
    apf_data.m->getDownward(el_apf, 2, down);
    apf::MeshEntity* face_apf = down[local_face_apf];

    face_weights[i] = apf_data.m->isGhost(face_apf) ? 0 : 1; 
  }
}


void getDofNums(ApfData& apf_data, const MeshEntityGroupSpec& vol_group,
                std::vector<apf::MeshEntity*>& elements,
                ArrayType<Index, 2>& nodenums)
{
  std::vector<int> nodenums_tmp;
  int el_idx = 0;
  for( auto e : elements )
  {
    getDofNums(apf_data, e, nodenums_tmp);
    for (unsigned int j=0; j < nodenums_tmp.size(); ++j)
      nodenums[el_idx][j] = nodenums_tmp[j];

    el_idx += 1;
  }
}

void getDofNums(ApfData& apf_data, apf::MeshEntity* e, std::vector<int>& node_nums)
{
  //node_nums.resize(0);
  node_nums.resize(apf_data.m_ref_el_sol->getNumNodesTotal());
  apf::Downward down;
  //std::vector<int> order(node_nums.size());
  for (int dim=0; dim <= 3; ++dim)
  {
    int ndown = apf_data.m->getDownward(e, dim, down);
    apf::Mesh::Type etype = apf_data.m->getType(down[0]);
    //apf::EntityShape* eshape = apf_data.sol_shape->getEntityShape(etype);
    int nnodes_dim = apf_data.sol_shape->countNodesOn(etype);
    for (int i=0; i < ndown; ++i)
    {
      int entity_ref_el = apf_data.m_ref_el_sol->getREEntityIndex(dim, i);
      //eshape->alignSharedNodes(apf_data.m, e, down[i], order.data());
      for (int j=0; j < nnodes_dim; ++j)
      {
        //TODO: this isn't right for 3rd order and above: the jth node on the apf entity
        //      might not be the jth node on the RE entity, because orientations might not
        //      be the same
        int apf_canonical_node = j; // order[j];
        int node = apf_data.m_ref_el_sol->getNodeIndex(dim, entity_ref_el, apf_canonical_node);
        assert(node >= 0 && node < node_nums.size());
        node_nums[node] = apf::getNumber(apf_data.dof_nums, down[i], j, 0);
      }
    }
  }
}


void getCoords(ApfData& apf_data, const MeshEntityGroupSpec& vol_group,
               std::vector<apf::MeshEntity*>& elements,
               ArrayType<Real, 3>& coords)
{
  apf::Downward down;
  apf::FieldShape* fshape = apf_data.m->getShape();
  apf::Vector3 point;
  int el_idx = 0;
  for( auto e : elements )
  {
    for (int dim=0; dim <= 3; ++dim)
    {
      int ndown = apf_data.m->getDownward(e, dim, down);

      int nnodes_dim = fshape->countNodesOn(apf_data.m->getType(down[0]));
      for (int i=0; i < ndown; ++i)
        for (int j=0; j < nnodes_dim; ++j)
        {
          apf_data.m->getPoint(down[i], j, point);

          int entity_ref_el = apf_data.m_ref_el_coord->getREEntityIndex(dim, i);
          int node = apf_data.m_ref_el_coord->getNodeIndex(dim, entity_ref_el, j);

          coords[el_idx][node][0] = point.x();
          coords[el_idx][node][1] = point.y();
          coords[el_idx][node][2] = point.z();
        }
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

      apf_data.m->end(it);
    }
}

/*
ArrayType<LocalIndex, 2> getFaceNodeMap(ApfData apf_data, apf::FieldShape* fshape)
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

      int nnodes_dim = fshape->countNodesOn(apf_data.m->getType(face_entities[0]));
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
*/
/*
const ArrayType<LocalIndex, 2>& getFaceTensorProductMap(const int degree)
{
  // this has to be consistent with getFaceNodemap
  // this only works for coordinate field
  static ArrayType<LocalIndex, 2> degree1(boost::extents[2][2]);
  static ArrayType<LocalIndex, 2> degree2(boost::extents[3][3]);

  degree1[0][0] = 0;
  degree1[1][0] = 1;
  degree1[1][1] = 2;
  degree1[0][1] = 3;

  degree2[0][0] = 0;
  degree2[2][0] = 1;
  degree2[2][2] = 2;
  degree2[0][2] = 3;
  degree2[1][0] = 4;
  degree2[2][1] = 5;
  degree2[1][2] = 6;
  degree2[0][1] = 7;
  degree2[1][1] = 8;

  if (degree == 1)
    return degree1;
  else if (degree == 2)
    return degree2;
  else
    throw std::invalid_argument("unsupported degree");
}
*/

} // namespace
