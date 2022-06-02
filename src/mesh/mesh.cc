#include "mesh/mesh.h"
#include "memory.h"
#include "mesh/apfMDSField.h"
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
  m_volume_spec(volume_group_spec),
  m_bc_spec(bc_spec),
  m_all_face_spec(bc_spec),
  m_ref_el_coord(reference_element::getLagrangeHexReferenceElement(coord_degree)),
  m_ref_el_sol(reference_element::getLagrangeHexReferenceElement(solution_degree)),
  m_apf_data(m, m_ref_el_sol, m_ref_el_coord),
  m_tensor_product_coord_map(m_ref_el_coord),
  m_tensor_product_sol_map(m_ref_el_sol),
  m_field_data_manager(m_apf_data.m, m_apf_data.sol_shape, m_apf_data.dof_nums, m_apf_data.is_dirichlet)
{
  assert(coord_degree == 1);
  m_dof_numbering.sol_degree   = solution_degree;
  m_dof_numbering.coord_degree = coord_degree;
  setVolumeIndices();
  setSurfaceIndices(other_surface_spec);
  setApfData();
  createVolumeGroups();
  createFaceGroups();
}


void MeshCG::setVolumeIndices()
{
  for (SInt i=0; i < m_volume_spec.size(); ++i)
    m_volume_spec[i].setIdx(i);
}

void MeshCG::setSurfaceIndices(const std::vector<MeshEntityGroupSpec>& other_surface_spec)
{
  // assign indices to surfaces, also copies the surfaces from bc_spec
  // into all_faces_spec (which initially contains only non-bc surfaces)

  for (auto& surf_spec : other_surface_spec)
  {
    m_all_face_spec.push_back(surf_spec);
    if (surf_spec.getIsDirichlet())
      throw std::runtime_error("non-BC surfaces cannot be marked as Dirichlet");
  }

  for (SInt i=0; i < m_bc_spec.size(); ++i)
    m_bc_spec[i].setIdx(i);

  for (SInt i=0; i < m_all_face_spec.size(); ++i)
    m_all_face_spec[i].setIdx(i);

  /*
  for (SInt i=0; i < m_all_face_spec.size(); ++i)
    m_all_face_spec[i].setIdx(i + m_bc_spec.size());

  for (SInt j=0; j < m_bc_spec.size(); ++j)
  {
    m_bc_spec[j].setIdx(j);
    m_all_face_spec.push_back(m_bc_spec[j]);
  }
  */
}




void MeshCG::setApfData()
{
  //TODO: make this the constructor of ApfData
  //m_apf_data.sol_shape = Mesh::getLagrange(m_dof_numbering.sol_degree);
  //m_apf_data.coord_shape = Mesh::getLagrange(m_dof_numbering.coord_degree);

  m_dof_numbering.nodes_per_element =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::HEX);
  m_dof_numbering.nodes_per_face =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::QUAD);
  m_dof_numbering.coord_nodes_per_element =
    apf::countElementNodes(m_apf_data.coord_shape, apf::Mesh::HEX);

  setVolumeGroupNumbering(m_apf_data.m, m_volume_spec, m_apf_data.vol_groups);

  //TODO: move dof numbering to here?
}

void MeshCG::createVolumeGroups()
{
  setDirichletNodes(m_apf_data, m_bc_spec);
  AdjacencyNumberer reorderer_local(m_apf_data.m, m_apf_data.dof_nums,
                                    m_apf_data.el_nums, m_apf_data.is_dirichlet);

  AdjacencyNumberer reorderer_global(m_apf_data.m, m_apf_data.global_dof_nums,
                                     nullptr, m_apf_data.is_dirichlet, true);
  reorderer_local.reorder();
  reorderer_global.reorder();
  int num_dofs_local = reorderer_local.getNumDofs();
  int num_dofs_owned = reorderer_global.getNumDofs();

  m_dof_numbering.num_dofs_local = num_dofs_local;
  m_dof_numbering.num_dofs_owned = num_dofs_owned;
  m_dof_numbering.num_dofs_total = reorderer_local.getNumTotalDofs();
  m_dof_numbering.local_owned_to_global_dof_offset = reorderer_global.getLocalOwnedToGlobalOffset();

  m_elements = reorderer_local.getElements();
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
    std::vector<int_least8_t> element_weights;

    getGroupElements(m_apf_data, m_volume_spec[i], elements_group);
    getElementWeights(m_apf_data, elements_group, element_weights);
    getDofNums(m_apf_data, m_volume_spec[i], elements_group, dof_nums);
    getCoords(m_apf_data,  m_volume_spec[i], elements_group, coords);

    m_vol_group.emplace_back(VolumeGroup(m_volume_spec[i].getIdx(), dof_nums, coords,
      m_tensor_product_coord_map, m_tensor_product_sol_map, 
      m_ref_el_coord, m_ref_el_sol, elements_group, element_weights));

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
  for (auto& surf : m_all_face_spec)
  {
    std::cout << "surface " << surf.getIdx() << " has model entities ";
    for (auto me : surf.getModelEntities())
      std::cout << me << ", ";
    std::cout << std::endl;

    bool is_boundary_surface = surf.getIdx() < m_bc_spec.size();
    std::vector<FaceSpec> faces;
    getGroupFaces(m_apf_data, surf, m_ref_el_coord, m_volume_spec, m_elnums_global_to_local, faces);

    std::vector<int_least8_t> face_weights;
    getFaceWeights(m_apf_data, faces, m_ref_el_coord, m_elements, face_weights);

    // get dofs
    auto nfaces = faces.size();
    std::cout << "nfaces = " << nfaces << std::endl;
    int num_nodes_per_face = m_ref_el_sol->getNumNodes(2);
    ArrayType<Index, 2> nodenums(boost::extents[nfaces][num_nodes_per_face]);
    auto& nodemap_sol = m_ref_el_sol->getFaceNodes();
    for (SInt i=0; i < nfaces; ++i)
    {
      FaceSpec& face_i = faces[i];
      int vol_group = apf::getNumber(m_apf_data.vol_groups, 
                                     m_elements[face_i.el], 0, 0);
      auto& parent_group = m_vol_group[vol_group];

      //TODO: consider storing VolumeGroup dofs in offset array (so
      //they can be indexed by the global element number rather than
      //the local one)
      int group_elnum = m_elnums_global_to_local[face_i.el];
      for (int node=0; node <  num_nodes_per_face; ++node)
        nodenums[i][node] = parent_group.nodenums[group_elnum][nodemap_sol[face_i.face][node]];
    }

    m_all_faces.emplace_back(surf.getIdx(), m_ref_el_coord, m_ref_el_sol,
                          m_tensor_product_coord_map, m_tensor_product_sol_map,
                          /*tp_nodemap,*/ faces, nodenums, face_weights, surf.getIsDirichlet(), is_boundary_surface);
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


// for non-owned dofs, gives both the global dof number and the corresponding local dof number
void MeshCG::getGhostDofInfo(std::vector<DofInt>& global_dofs, std::vector<DofInt>& local_dofs)
{
  int ncomp = apf::countComponents(m_apf_data.dof_nums);
  std::vector<std::pair<DofInt, DofInt>> data;

  for (int dim=0; dim <= 3; ++dim)
  {
    auto it = m_apf_data.m->begin(dim);
    apf::MeshEntity* e;
    while ( (e = m_apf_data.m->iterate(it)))
    {
      if (!m_apf_data.m->isOwned(e))
      {
        apf::Mesh::Type etype = m_apf_data.m->getType(e);
        int nnodes_dim = m_apf_data.sol_shape->countNodesOn(etype);
        for (int i=0; i < nnodes_dim; ++i)
          for (int c=0; c < ncomp; ++c)
          {
            int local_dof_num = apf::getNumber(m_apf_data.dof_nums, e, i, c);
            if (isDofActive(local_dof_num))
            {
              int global_dof_num = apf::getNumber(m_apf_data.global_dof_nums, e, i, c);
              data.emplace_back(local_dof_num, global_dof_num);

            }
          }
      }
    }
    m_apf_data.m->end(it);
  }

  //sort both vectors by global_dofs (which should correspond to an adjacency-based ordering), and should be
  // cache friendly
  auto compare_global = [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return std::less<int>()(p1.second, p2.second); };
  std::sort(data.begin(), data.end(), compare_global);
  global_dofs.resize(data.size());
  local_dofs.resize(data.size());
  for (size_t i=0; i < data.size(); ++i)
  {
    global_dofs[i] = data[i].second;
    local_dofs[i]  = data[i].first;
  }
}

void MeshCG::getOwnedLocalDofInfo(std::vector<DofInt>& owned_local_dofs)
{
  owned_local_dofs.resize(getNumOwnedDofs());
  int global_to_local_offset = m_dof_numbering.local_owned_to_global_dof_offset;
  int ncomp = apf::countComponents(m_apf_data.dof_nums);

  for (int dim=0; dim <= 3; ++dim)
  {
    auto it = m_apf_data.m->begin(dim);
    apf::MeshEntity* e;
    while ( (e = m_apf_data.m->iterate(it)))
    {
      if (m_apf_data.m->isOwned(e))
      {
        apf::Mesh::Type etype = m_apf_data.m->getType(e);
        int nnodes_dim = m_apf_data.sol_shape->countNodesOn(etype);
        for (int i=0; i < nnodes_dim; ++i)
          for (int c=0; c < ncomp; ++c)
          {
            int local_dof_num = apf::getNumber(m_apf_data.dof_nums, e, i, c);
            if (isDofActive(local_dof_num))
            {
              int local_owned_dof = apf::getNumber(m_apf_data.global_dof_nums, e, i, c) - global_to_local_offset;
              owned_local_dofs[local_owned_dof] = local_dof_num;
            }
          }
      }
    }

    m_apf_data.m->end(it);
  }
}


void MeshCG::getLocalToGlobalDofs(std::vector<DofInt>& local_to_global_dofs)
{
  local_to_global_dofs.resize(getNumDofs());
  int ncomp = apf::countComponents(m_apf_data.dof_nums);

  for (int dim=0; dim <= 3; ++dim)
  {
    auto it = m_apf_data.m->begin(dim);
    apf::MeshEntity* e;
    while ( (e = m_apf_data.m->iterate(it)))
    {
      if (m_apf_data.m->isOwned(e))
      {
        apf::Mesh::Type etype = m_apf_data.m->getType(e);
        int nnodes_dim = m_apf_data.sol_shape->countNodesOn(etype);
        for (int i=0; i < nnodes_dim; ++i)
          for (int c=0; c < ncomp; ++c)
          {
            int local_dof_num = apf::getNumber(m_apf_data.dof_nums, e, i, c);
            if (isDofActive(local_dof_num))
            {
              int local_dof = apf::getNumber(m_apf_data.dof_nums, e, i, c);
              local_to_global_dofs[local_dof] = apf::getNumber(m_apf_data.global_dof_nums, e, i, c);
            }
          }
      }
    }

    m_apf_data.m->end(it);
  }
}

std::shared_ptr<MeshCG> createMeshCG(apf::Mesh2* m,
                                     std::vector<MeshEntityGroupSpec> volume_group_spec,
                                     std::vector<MeshEntityGroupSpec> bc_spec,
                                     std::vector<MeshEntityGroupSpec> other_surface_spec,
                                     const int solution_degree, const int coord_degree)
{
  createGhostLayer(m);
  return std::shared_ptr<MeshCG>(new MeshCG(m, volume_group_spec, bc_spec, other_surface_spec, solution_degree, coord_degree));
}





} // namespace
