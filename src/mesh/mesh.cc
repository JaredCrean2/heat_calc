#include "mesh/mesh.h"
#include "mesh/mesh_helper.h"
#include "mesh/dof_numbering.h"
#include "mesh/apfShapeHex.h"
#include "PCU.h"
#include "mpi.h"

#include "apfShape.h"

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
  m_tensor_product_coord_map(getReferenceElement(apf::Mesh::HEX, coord_degree)),
  m_tensor_product_sol_map(getReferenceElement(apf::Mesh::HEX, solution_degree))
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

  m_dof_numbering.nodes_per_element =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::HEX);
  m_dof_numbering.nodes_per_face =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::QUAD);
  m_dof_numbering.coord_nodes_per_element =
    apf::countElementNodes(m_apf_data.coord_shape, apf::Mesh::HEX);

  m_apf_data.dof_nums = apf::createNumbering(m_apf_data.m, "dof_nums",
                                          m_apf_data.sol_shape, 1);
  m_apf_data.el_nums = apf::createNumbering(m_apf_data.m, "el_nums",
                                        apf::getConstant(3), 1);
  m_apf_data.is_dirichlet = apf::createNumbering(m_apf_data.m, "is_dirichlet",
                                              m_apf_data.sol_shape, 1);
  m_apf_data.vol_groups = apf::createNumbering(m_apf_data.m, "vol_group",
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

  ReferenceElement* ref_el_coord = getReferenceElement(apf::Mesh::HEX,
                                            m_dof_numbering.coord_degree);
  ReferenceElement* ref_el_sol   = getReferenceElement(apf::Mesh::HEX,
                                            m_dof_numbering.sol_degree);

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
      ref_el_coord, ref_el_sol, elements_group));

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
  ArrayType<LocalIndex, 2> nodemap_coord = getFaceNodeMap(m_apf_data, m_apf_data.coord_shape);
  ArrayType<LocalIndex, 2> nodemap_sol = getFaceNodeMap(m_apf_data, m_apf_data.sol_shape);
  const auto& tp_nodemap = getFaceTensorProductMap(m_dof_numbering.coord_degree);

  ReferenceElement* ref_el_coord = getReferenceElement(apf::Mesh::HEX,
                                            m_dof_numbering.coord_degree);
  ReferenceElement* ref_el_sol   = getReferenceElement(apf::Mesh::HEX,
                                            m_dof_numbering.sol_degree);

  for (auto& surf : m_all_face_spec)
  {
    m_all_faces.emplace_back(surf.getIdx(), ref_el_coord, ref_el_sol,
                             m_tensor_product_coord_map, m_tensor_product_sol_map,
                             nodemap_coord, nodemap_sol,
                             tp_nodemap, surf.getIsDirichlet());
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
    int num_nodes_per_face = nodemap_sol.shape()[1];
    face_group.nodenums.resize(boost::extents[nfaces][num_nodes_per_face]);
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

} // namespace
