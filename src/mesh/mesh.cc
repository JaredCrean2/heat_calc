#include "mesh/mesh.h"
#include "mesh/mesh_helper.h"
#include "mesh/dof_numbering.h"
#include "mesh/apfShapeHex.h"
#include "PCU.h"
#include "mesh/reference_element_apf.h"
#include "mesh/reference_element_interface.h"
#include "mpi.h"

#include "apfShape.h"
#include <apfNumbering.h>

namespace Mesh {

bool _initialized = initialize();

bool initialize()
{
  if (!PCU_Comm_Initialized())
  {
    MPI_Init(0, NULL);
    PCU_Comm_Init();
  }

  return true;
}

MeshCG::MeshCG(apf::Mesh2* m,
        std::vector<MeshEntityGroupSpec> volume_group_spec,
        std::vector<MeshEntityGroupSpec> bc_spec,
        std::vector<MeshEntityGroupSpec> other_surface_spec,
        const int solution_degree, const int coord_degree) :
  m_apf_data(m),
  m_volume_spec(volume_group_spec),
  m_bc_spec(bc_spec),
  m_all_face_spec(other_surface_spec),
  m_ref_el_coord(reference_element::getLagrangeHexReferenceElement(coord_degree)),
  m_ref_el_sol(reference_element::getLagrangeHexReferenceElement(solution_degree)),
  m_tensor_product_coord_map(m_ref_el_coord),
  m_tensor_product_sol_map(m_ref_el_sol)
{
  assert(coord_degree == 1);
  m_dof_numbering.sol_degree   = solution_degree;
  m_dof_numbering.coord_degree = coord_degree;
  setVolumeIndices();
  setSurfaceIndices();
  setApfData();
  createVolumeGroups();
  createFaceGroups();
}


void MeshCG::setVolumeIndices()
{
  for (SInt i=0; i < m_volume_spec.size(); ++i)
    m_volume_spec[i].setIdx(i);
}

void MeshCG::setSurfaceIndices()
{
  // assign indices to surfaces, also copies the surfaces from bc_spec
  // into all_faces_spec (which initially contains only non-bc surfaces)
  for (SInt i=0; i < m_all_face_spec.size(); ++i)
    m_all_face_spec[i].setIdx(i + m_bc_spec.size());

  for (SInt j=0; j < m_bc_spec.size(); ++j)
  {
    m_bc_spec[j].setIdx(j);
    m_all_face_spec.push_back(m_bc_spec[j]);
  }
}


void MeshCG::setApfData()
{
  m_apf_data.sol_shape = Mesh::getLagrange(m_dof_numbering.sol_degree);
  m_apf_data.coord_shape = Mesh::getLagrange(m_dof_numbering.coord_degree);

  auto sol_shape   = apf::getHexFieldShape(m_ref_el_sol);
  auto coord_shape = apf::getHexFieldShape(m_ref_el_coord);
  m_apf_data.m_sol_shape    = sol_shape;
  m_apf_data.sol_shape      = sol_shape.get();
  m_apf_data.m_coord_shape  = coord_shape;
  m_apf_data.coord_shape    = coord_shape.get();
  m_apf_data.m_ref_el_coord = m_ref_el_coord;
  m_apf_data.m_ref_el_sol   = m_ref_el_sol;

  m_dof_numbering.nodes_per_element =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::HEX);
  m_dof_numbering.nodes_per_face =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::QUAD);
  m_dof_numbering.coord_nodes_per_element =
    apf::countElementNodes(m_apf_data.coord_shape, apf::Mesh::HEX);

  m_apf_data.dof_nums = apf::createNumberingMDS(m_apf_data.m, "dof_nums",
                                          m_apf_data.sol_shape, 1);
  m_apf_data.el_nums = apf::createNumberingMDS(m_apf_data.m, "el_nums",
                                        apf::getConstant(3), 1);
  m_apf_data.is_dirichlet = apf::createNumberingMDS(m_apf_data.m, "is_dirichlet",
                                              m_apf_data.sol_shape, 1);
  m_apf_data.vol_groups = apf::createNumberingMDS(m_apf_data.m, "vol_group",  //TODO: is this needed?
                                               apf::getConstant(3), 1);

  setVolumeGroupNumbering(m_apf_data.m, m_volume_spec, m_apf_data.vol_groups);

}

void MeshCG::createVolumeGroups()
{
  setDirichletNodes(m_apf_data, m_bc_spec);
  AdjacencyNumberer reorderer(m_apf_data.m, m_apf_data.dof_nums,
                              m_apf_data.el_nums, m_apf_data.is_dirichlet);
  reorderer.reorder();
  int num_dofs = reorderer.getNumDofs();

  m_dof_numbering.num_dofs = num_dofs;
  m_dof_numbering.num_dofs_total = reorderer.getNumTotalDofs();
  //apf::synchronize(m_apf_data.dof_nums);

  m_elements = reorderer.getElements();
  m_elnums_global_to_local.resize(m_elements.size());
  m_elnums_local_to_global.resize(m_volume_spec.size());

  for (SInt i=0; i < m_volume_spec.size(); ++i)
  {
    //TODO: avoid Boost array deep copies
    int numel_i = countNumEls(m_apf_data, m_volume_spec[i]);
    ArrayType<Index, 2> dof_nums(boost::extents[numel_i][m_dof_numbering.nodes_per_element]);
    ArrayType<Real, 3> coords(boost::extents[numel_i][m_dof_numbering.coord_nodes_per_element][3]);
    //TODO: get coords on coordinate field
    std::vector<apf::MeshEntity*> elements_group;

    getGroupElements(m_apf_data, m_volume_spec[i], elements_group);
    getDofNums(m_apf_data, m_volume_spec[i], elements_group, dof_nums);
    getCoords(m_apf_data,  m_volume_spec[i], elements_group, coords);

    m_vol_group.emplace_back(VolumeGroup(m_volume_spec[i].getIdx(), dof_nums, coords,
      m_tensor_product_coord_map, m_tensor_product_sol_map, 
      m_ref_el_coord, m_ref_el_sol, elements_group));

    m_elnums_local_to_global[i].resize(elements_group.size());
    for (SInt j=0; j < elements_group.size(); ++j)
    {
      Index elnum_global = apf::getNumber(m_apf_data.el_nums, elements_group[j], 0, 0);
      m_elnums_global_to_local[elnum_global] = j;
      m_elnums_local_to_global[i][j] = elnum_global;
    }
  }

}

void MeshCG::createFaceGroups()
{

  apf::Downward down;
  //ArrayType<LocalIndex, 2> nodemap_coord = getFaceNodeMap(m_apf_data, m_apf_data.coord_shape);
  //ArrayType<LocalIndex, 2> nodemap_sol = getFaceNodeMap(m_apf_data, m_apf_data.sol_shape);
  //const auto& tp_nodemap = getFaceTensorProductMap(m_dof_numbering.coord_degree);
/*
  ReferenceElement* ref_el_coord = getReferenceElement(apf::Mesh::HEX,
                                            m_dof_numbering.coord_degree);
  ReferenceElement* ref_el_sol   = getReferenceElement(apf::Mesh::HEX,
                                            m_dof_numbering.sol_degree);
*/
  for (auto& surf : m_all_face_spec)
  {
    m_all_faces.emplace_back(surf.getIdx(), m_ref_el_coord, m_ref_el_sol,
                             m_tensor_product_coord_map, m_tensor_product_sol_map,
                             /*tp_nodemap,*/ surf.getIsDirichlet());
    auto& face_group = m_all_faces.back();

    //TODO: consider doing adjacency based search (starting with min
    // element number of mirror dof numbering algorithm)
    apf::MeshIterator* it = m_apf_data.m->begin(2);
    apf::MeshEntity* e;
    while ( (e = m_apf_data.m->iterate(it)) )
    {
      auto me = getMESpec(m_apf_data.m, e);
      if (surf.hasModelEntity(me))
      {
        auto me_parent = surf.getParentEntity(me);
        for (int i=0; i < m_apf_data.m->countUpward(e); ++i)
        {
          apf::MeshEntity* e_i = m_apf_data.m->getUpward(e, i);
          auto me_i = getMESpec(m_apf_data.m, e_i);
          if (me_i == me_parent)
          {
            // get local face number
            int nfaces_i = m_apf_data.m->getDownward(e_i, 2, down);
            int localidx = -1;
            for (int j=0; j < nfaces_i; ++j)
              if (down[j] == e)
              {
                localidx = j;
                break;
              }

            assert(localidx != -1);
            localidx = m_ref_el_coord->getREEntityIndex(2, localidx);

            int elnum = apf::getNumber(m_apf_data.el_nums, e_i, 0, 0);
            int elnum_local = m_elnums_global_to_local[elnum];
            int vol_group_idx = -1;
            // TODO: could cache this, because it is a function of me_parent
            for (const auto& vol_group : m_volume_spec)
              if (vol_group.hasModelEntity(me_parent))
                vol_group_idx = vol_group.getIdx();

            assert(vol_group_idx != -1);

            face_group.faces.push_back(FaceSpec(elnum, elnum_local, localidx, vol_group_idx));
            break;
          }
        }
      }
    }
    m_apf_data.m->end(it);


    // get dofs
    auto nfaces = face_group.faces.size();
    int num_nodes_per_face = m_ref_el_sol->getNumNodes(2);
    face_group.nodenums.resize(boost::extents[nfaces][num_nodes_per_face]);
    auto& nodemap_sol = m_ref_el_sol->getFaceNodes();
    //face_group.nodenums = ArrayType<Index, 2>(boost::extents[nfaces][num_nodes_per_face]);

    for (SInt i=0; i < nfaces; ++i)
    {
      FaceSpec& face_i = face_group.faces[i];
      int vol_group = apf::getNumber(m_apf_data.vol_groups, 
                                     m_elements[face_i.el], 0, 0);
      auto& parent_group = m_vol_group[vol_group];

      //TODO: consider storing VolumeGroup dofs in offset array (so
      //they can be indexed by the global element number rather than
      //the local one)
      int group_elnum = m_elnums_global_to_local[face_i.el];
      for (int node=0; node <  num_nodes_per_face; ++node)
        face_group.nodenums[i][node] =
          parent_group.nodenums[group_elnum][nodemap_sol[face_i.face][node]];
    }
  }
}


void MeshCG::getElementDofs(const VolumeGroup& vol_group, int el_idx, std::vector<Index>& nodenums)
{
  int elnum_global = m_elnums_local_to_global[vol_group.getIdx()].at(el_idx);
  getDofNums(m_apf_data, m_elements[elnum_global], nodenums);
}

/*
void MeshCG::getDofConnectivity(const VolumeGroup& vol_group, const Index el_idx, const Index node, std::vector<DofInt>& dofs)
{

  auto node_location = m_apf_data.m_ref_el_sol->getNodeLocation(node);
  int entity_apf = m_apf_data.m_ref_el_sol->getApfEntityIndex(node_location.dim, node_location.entity_index);
  apf::MeshEntity* el = vol_group.m_elements[el_idx];
  apf::Downward down;
#ifdef NDEBUG
  m_apf_data.m->getDownward(el, node_location.dim, down);
#else
  int ndown = m_apf_data.m->getDownward(el, node_location.dim, down);
  assert(entity_apf < ndown);
#endif
  apf::MeshEntity* entity = down[entity_apf];

  apf::Adjacent other_els;
  m_apf_data.m->getAdjacent(entity, 3, other_els);

  dofs.resize(0);
  std::vector<DofInt> dofs_el;
  for (auto& other_el : other_els)
  {
    getDofNums(m_apf_data, other_el, dofs_el);
    for (auto& dof : dofs_el)
      dofs.push_back(dof);
  }

  std::sort(dofs.begin(), dofs.end());
  auto it = std::unique(dofs.begin(), dofs.end());
  dofs.erase(it, dofs.end());
}
*/

void MeshCG::getDofConnectivity(const VolumeGroup& vol_group, const Index el_idx, const Index node, std::vector<DofInt>& dofs)
{
  dofs.resize(0);
  auto node_location = m_apf_data.m_ref_el_sol->getNodeLocation(node);
  int entity_apf = m_apf_data.m_ref_el_sol->getApfEntityIndex(node_location.dim, node_location.entity_index);
  apf::MeshEntity* el = vol_group.m_elements[el_idx];
  apf::Downward down;
#ifdef NDEBUG
  m_apf_data.m->getDownward(el, node_location.dim, down);
#else
  int ndown = m_apf_data.m->getDownward(el, node_location.dim, down);
  assert(entity_apf < ndown);
#endif
  apf::MeshEntity* entity = down[entity_apf];

  apf::Adjacent other_els;
  m_apf_data.m->getAdjacent(entity, 3, other_els);

  std::array<std::vector<apf::MeshEntity*>, 4> unique_entities;
  for (int dim=0; dim <= 3; ++dim)
    unique_entities[dim].reserve(12*8);

  for (auto& other_el : other_els)
    for (int dim=0; dim <= 3; ++dim)
    {
      if (!m_apf_data.sol_shape->hasNodesIn(dim))
        break;

      int ndown = m_apf_data.m->getDownward(other_el, dim, down);
      for (int i=0; i < ndown; ++i)
        unique_entities[dim].push_back(down[i]);
    }

  for (int dim=0; dim <= 3; ++dim)
  {
    if (!m_apf_data.sol_shape->hasNodesIn(dim))
      break;

    auto& entities_dim = unique_entities[dim];
    std::sort(entities_dim.begin(), entities_dim.end());
    auto it_end = std::unique(entities_dim.begin(), entities_dim.end());

    //for (auto& e : unique_entities[dim])
    for (auto it=entities_dim.begin(); it != it_end; ++it)
      for (int j=0; j < m_ref_el_sol->getNumNodes(dim); ++j)
        dofs.push_back(apf::getNumber(m_apf_data.dof_nums, *it, j, 0));
  }
}



} // namespace
