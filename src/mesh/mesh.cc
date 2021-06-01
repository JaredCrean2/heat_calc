#include "mesh/mesh.h"
#include "mesh/mesh_helper.h"
#include "apfShape.h"

namespace Mesh {

MeshCG::MeshCG(apf::Mesh* m,
        std::vector<MeshEntityGroupSpec> volume_group_spec,
        std::vector<MeshEntityGroupSpec> bc_spec,
        std::vector<MeshEntityGroupSpec> other_surface_spec,
        int solution_degree, int coord_degree) :
  m_apf_data(m),
  m_volume_spec(volume_group_spec),
  m_bc_spec(bc_spec),
  m_all_face_spec(other_surface_spec)
{
  assert(coord_degree == 1);
  m_dof_numbering.sol_degree   = solution_degree;
  m_dof_numbering.coord_degree = coord_degree;
  setSurfaceIndices();
  setApfData();
  createVolumeGroups();
  createFaceGroups();
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
  m_apf_data.sol_shape = apf::getLagrange(m_dof_numbering.sol_degree);
  m_dof_numbering.nodes_per_element =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::HEX);
  m_dof_numbering.nodes_per_face =
    apf::countElementNodes(m_apf_data.sol_shape, apf::Mesh::QUAD);

  m_apf_data.dof_nums = apf::createNumbering(m_apf_data.m, "dof_nums",
                                          m_apf_data.sol_shape, 1);
  m_apf_data.el_nums = apf::createNumbering(m_apf_data.m, "el_nums",
                                        apf::getConstant(3), 1);
  m_apf_data.is_dirichlet = apf::createNumbering(m_apf_data.m, "is_dirichlet",
                                              m_apf_data.sol_shape, 1);
}

void MeshCG::createVolumeGroups()
{
  setDirichletNodes(m_apf_data, m_bc_spec);

  int elnum_start = 0, elnum_end = 0, curr_dof = 0;
  std::vector<std::pair<int, int>> el_bounds;
  initializeElnums(m_apf_data);
  for (SInt i=0; i < m_volume_spec.size(); ++i)
  {
    setDofNums(m_apf_data.m, m_volume_spec[i],
               elnum_start, curr_dof,
               m_apf_data.is_dirichlet,
               m_apf_data.dof_nums, m_apf_data.el_nums);
    el_bounds.push_back(std::make_pair(elnum_start, elnum_end-1));
    elnum_start = elnum_end;
  }

  m_dof_numbering.num_dofs = curr_dof - 1;
  setDirichletDofs(m_apf_data, curr_dof);
  //apf::synchronize(m_apf_data.dof_nums);
  m_dof_numbering.num_dofs_total = curr_dof - 1;

  for (SInt i=0; i < m_volume_spec.size(); ++i)
  {
    auto bounds = el_bounds[i];
    //TODO: avoid Boost array deep copies
    ArrayType<Index, 2> dof_nums(boost::extents[bounds.second - bounds.first + 1][m_dof_numbering.nodes_per_element]);
    getDofNums(m_apf_data, bounds.first, bounds.second, dof_nums);
    m_vol_group.push_back(VolumeGroup(bounds.first, bounds.second, dof_nums));
  }
}

void MeshCG::createFaceGroups()
{

  apf::Downward down;
  ArrayType<LocalIndex, 2> nodemap = getFaceNodeMap(m_apf_data);
  for (auto& surf : m_all_face_spec)
  {
    m_all_faces.emplace_back();
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

            int elnum = apf::getNumber(m_apf_data.el_nums, e, 0, 0);
            face_group.faces.push_back(FaceSpec(elnum, localidx));
          }
        }
      }
    }


    // get dofs
    auto nfaces = face_group.faces.size();
    int num_nodes_per_face = nodemap.shape()[1];
    face_group.nodenums = ArrayType<Index, 2>(boost::extents[nfaces][num_nodes_per_face]);

    for (SInt i=0; i < nfaces; ++i)
    {
      FaceSpec& face_i = face_group.faces[i];

      // get VolumeGroup
      VolumeGroup* parent_group= nullptr;
      for (auto& volgroup : m_vol_group)
        if (face_i.el >= volgroup.el_start &&
            face_i.el <= volgroup.el_end)
          parent_group = &volgroup;

      assert(parent_group != nullptr);
      //TODO: consider storing VolumeGroup dofs in offset array (so
      //they can be indexed by the global element number rather than
      //the local one)
      int group_elnum = face_i.el - parent_group->el_start; 
      for (int node=0; node <  num_nodes_per_face; ++node)
        face_group.nodenums[i][node] =
          parent_group->nodenums[group_elnum][nodemap[face_i.face][node]];
    }
  }
}

} // namespace
