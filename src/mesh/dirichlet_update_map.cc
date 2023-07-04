#include "mesh/dirichlet_update_map.h"
#include "mesh/mesh.h"
#include "parallel_exchange.h"

namespace Mesh {

DirichletUpdateMap::DirichletUpdateMap(MeshCG* mesh, MPI_Comm comm) :
  m_mesh(mesh),
  m_apf_data(mesh->m_apf_data),
  m_comm(comm),
  m_send_nodes(commSize(comm)),
  m_recv_counts(commSize(comm), 0),
  m_recv_nodes(commSize(comm)),
  m_recv_node_ptrs(commSize(comm))
{
  getLocalDirichletUpdateMap();
  auto min_dirichlet_surface = collectMinDirichletSurfaces();
  createSendAndRecvLists(min_dirichlet_surface);
}


std::shared_ptr<DirichletUpdateMap::ModelEntityField> DirichletUpdateMap::collectMinDirichletSurfaces()
{
  ParallelExchange<SurfaceModelEntity> exchanger(MPI_COMM_WORLD, 1100);
  auto min_dirichlet_surface = std::make_shared<ModelEntityField> (m_apf_data.m, 
                                  "min_dirichlet_surface", apf::getLagrange(1), 1);

  apf::Copies remotes;
  apf::MeshIterator* it = m_apf_data.m->begin(0);
  apf::MeshEntity* e;
  while ( (e = m_apf_data.m->iterate(it)) )
  {
    if (apf::getNumber(m_apf_data.is_dirichlet, e, 0, 0) && isSharedOrGhost(e))
    {
      SurfaceModelEntity surf = getMinDirichletSurface(e);
      (*min_dirichlet_surface)(e, 0, 0) = surf;

      remotes.clear();
      m_apf_data.m->getRemotes(e, remotes);
      m_apf_data.m->getGhosts(e, remotes);

      for (auto& p : remotes)
      {
        surf.vert = p.second;
        exchanger.getSendBuf(p.first).push_back(surf);
        exchanger.getRecvBuf(p.first).emplace_back();
      }
    }
  }
  m_apf_data.m->end(it);

  exchanger.startCommunication();

  auto unpacker = [&](int rank, const std::vector<SurfaceModelEntity>& buf)
  {
    for (const auto& remote_entity : buf)
    {
      apf::MeshEntity* vert = remote_entity.vert;
      SurfaceModelEntity& local_entity = (*min_dirichlet_surface)(vert, 0, 0);
      if  (remote_entity < local_entity)
      {
        local_entity.rank = rank;
        local_entity.model_entity_dim = remote_entity.model_entity_dim;
        local_entity.model_entity_tag = remote_entity.model_entity_tag;
      }
    }
  };

  exchanger.finishCommunication(unpacker);

  return min_dirichlet_surface;
}

DirichletUpdateMap::SurfaceModelEntity DirichletUpdateMap::getMinDirichletSurface(apf::MeshEntity* e)
{
  SurfaceModelEntity min_surface;
  min_surface.rank = commRank(m_comm);
  apf::Adjacent faces;
  m_apf_data.m->getAdjacent(e, 2, faces);

  for (size_t i=0; i < faces.getSize(); ++i)
  {
    apf::MeshEntity* face = faces[i];
    apf::ModelEntity* me = m_apf_data.m->toModel(face);
    int dim = m_apf_data.m->getModelType(me);
    int tag = m_apf_data.m->getModelTag(me);

    if (dim == 2)
    {
      if (isSurfaceDirichlet(ModelEntitySpec{dim, tag}))
      {
        if (tag < min_surface.model_entity_tag)
        {
          min_surface.model_entity_dim = dim;
          min_surface.model_entity_tag = tag;
        }
      }
    }
  }

  return min_surface;
}

bool DirichletUpdateMap::hasLocalDirichletSurface(apf::MeshEntity* e)
{
  SurfaceModelEntity min_surface;
  min_surface.rank = commRank(m_comm);
  apf::Adjacent faces;
  m_apf_data.m->getAdjacent(e, 2, faces);

  for (size_t i=0; i < faces.getSize(); ++i)
  {
    apf::MeshEntity* face = faces[i];
    apf::ModelEntity* me = m_apf_data.m->toModel(face);
    int dim = m_apf_data.m->getModelType(me);
    int tag = m_apf_data.m->getModelTag(me);

    if (dim == 2 && isSurfaceDirichlet(ModelEntitySpec{dim, tag}))
    {
      return true;
    }
  }

  return false;
}

bool DirichletUpdateMap::isSurfaceDirichlet(const ModelEntitySpec& me_spec)
{
  for (auto& spec : m_mesh->m_bc_spec)
    if (spec.getIsDirichlet() && spec.hasModelEntity(me_spec))
      return true;

  return false;
}

void DirichletUpdateMap::createSendAndRecvLists(std::shared_ptr<ModelEntityField> min_dirichlet_surface)
{
  // loop over entities
  // skip if not dirichlet or shared
  // if have local dirichlet source, send false to min_dirichlet_surface, otherwise true
  ParallelExchange<MeshEntityStatus> exchanger(m_comm, 1102);
  for (auto& recv_node_ptrs : m_recv_node_ptrs)
    recv_node_ptrs.push_back(0);

  int myrank = commRank(m_comm);
  apf::Copies remotes;
  apf::MeshIterator* it = m_apf_data.m->begin(0);
  apf::MeshEntity* e;
  while ( (e = m_apf_data.m->iterate(it)) ) 
  {
    std::cout << "e = " << e << std::endl;
    std::cout << std::boolalpha;
    std::cout << "is dirichlet = " << bool(apf::getNumber(m_apf_data.is_dirichlet, e, 0, 0)) << ", isShared = " << bool(m_apf_data.m->isShared(e)) << ", isGhost = " << bool(m_apf_data.m->isGhost(e)) << ", isGhosted = " << bool(m_apf_data.m->isGhosted(e)) << std::endl;
    apf::Copies ghosts;
    m_apf_data.m->getGhosts(e, ghosts);
    for (auto& p : ghosts)
      std::cout << "  ghost " << p.second << " on rank " << p.first << std::endl;
    if (apf::getNumber(m_apf_data.is_dirichlet, e, 0, 0) && isSharedOrGhost(e))
    {
      std::cout << "min entity = " << (*min_dirichlet_surface)(e, 0, 0).model_entity_dim << ", " << (*min_dirichlet_surface)(e, 0, 0).model_entity_tag  << " on rank " << (*min_dirichlet_surface)(e, 0, 0).rank << std::endl;
      remotes.clear();      
      bool is_local = (*min_dirichlet_surface)(e, 0, 0).rank == myrank;
      bool have_local_source = hasLocalDirichletSurface(e);
      m_apf_data.m->getRemotes(e, remotes);
      m_apf_data.m->getGhosts(e, remotes);

      if (is_local)
      {
        for (auto& p : remotes)
        {
          std::cout << "  expecting to receive from rank " << p.first << ", entity " << p.second << std::endl;
          exchanger.getRecvBuf(p.first).push_back({nullptr, false});
        }        
      } else
      {
        int sender_rank = (*min_dirichlet_surface)(e, 0, 0).rank;
        std::cout << "  sending to rank " << sender_rank << " entity " << remotes[sender_rank] << std::endl;
        exchanger.getSendBuf(sender_rank).push_back({remotes[sender_rank], have_local_source});
        if (!have_local_source)
        {
          m_recv_counts[sender_rank]++;  //TODO: is this needed?
          getArrayNodes(e, m_recv_nodes[sender_rank]);
          m_recv_node_ptrs[sender_rank].push_back(m_recv_nodes[sender_rank].size());
        }
      }
    }
  }
  m_apf_data.m->end(it);

  std::cout << "myrank = " << commRank(m_comm) << std::endl;
  for (int i=0; i < commSize(m_comm); ++i)
    std::cout << "sending " << exchanger.getSendBuf(i).size() << " values to rank " << i << std::endl;

  for (int i=0; i < commSize(m_comm); ++i)
    std::cout << "receiving " << exchanger.getRecvBuf(i).size() << " values to rank " << i << std::endl;

  exchanger.startCommunication();

  auto unpacker = [&](int rank, const std::vector<MeshEntityStatus>& buf)
  {
    for (const auto& entity_status : buf)
    {
      if (!entity_status.status)
      {
        ArrayNode node = getArrayNode(entity_status.e);
        m_send_nodes[rank].push_back(node);
      }
    }
  };

  exchanger.finishCommunication(unpacker);
}

/*
void DirichletUpdateMap::createSendAndRecvLists(std::shared_ptr<ModelEntityField> min_dirichlet_surface)
{
  ParallelExchange<apf::MeshEntity*> exchanger(m_comm, 1101);
  int myrank = commRank(m_comm);
  apf::Copies remotes;
  apf::MeshIterator* it = m_apf_data.m->begin(0);
  apf::MeshEntity* e;
  while ( (e = m_apf_data.m->iterate(it)) )
  {
    if (apf::getNumber(m_apf_data.is_dirichlet, e, 0, 0) && m_apf_data.m->isShared(e))
    {
      bool is_local = (*min_dirichlet_surface)(e, 0, 0).rank == myrank;
      if (is_local)
      {
        m_apf_data.m->getRemotes(e, remotes);
        ArrayNode node = getArrayNode(e);
        for (auto& p : remotes)
        {
          m_send_nodes[p.first].push_back(node);
          exchanger.getSendBuf(p.first).push_back(p.second);
        }
      } else
      {
        int sender_rank = (*min_dirichlet_surface)(e, 0, 0).rank;
        exchanger.getRecvBuf(sender_rank).push_back(nullptr);
        m_recv_counts[sender_rank]++;
      }
    }  
  }
  m_apf_data.m->end(it);

  unpackRecvEntities(exchanger);
}

void DirichletUpdateMap::unpackRecvEntities(ParallelExchange<apf::MeshEntity*>& exchanger)
{
  auto unpacker = [&](int rank, const std::vector<apf::MeshEntity*>& buf)
  {
    auto& recv_nodes = m_recv_nodes[rank];
    auto& recv_nodes_ptr = m_recv_node_ptrs[rank];

    recv_nodes_ptr.push_back(0);
    std::vector<ArrayNode> nodes;
    for (apf::MeshEntity* vert : buf)
    {
      getArrayNodes(vert, recv_nodes);
      int end_idx = recv_nodes.size();
      recv_nodes_ptr.push_back(end_idx);
    }
  };

  exchanger.finishCommunication(unpacker);
}
*/

void DirichletUpdateMap::getLocalDirichletUpdateMap()
{
  apf::Adjacent els;
  m_local_dest_node_ptrs.push_back(0);
  for (auto& face_group : m_mesh->m_all_faces)
  {
    if (face_group.getIsDirichlet())
    {
      for (int i=0; i < face_group.getNumFaces(); ++i)
      {
        auto& face_spec = face_group.faces[i];
        apf::MeshEntity* el = m_mesh->m_elements[face_spec.el];

        apf::Downward faces, verts;
        m_apf_data.m->getDownward(el, 2, faces);
        int nverts = m_apf_data.m->getDownward(faces[m_mesh->m_ref_el_coord->getApfEntity(2, face_spec.face)], 0, verts);

        apf::ModelEntity* me = m_apf_data.m->toModel(el);

        for (int j=0; j < nverts; ++j)
        {
          Index start_idx = m_local_dest_nodes.size();
          Index end_idx = start_idx;
          m_apf_data.m->getAdjacent(verts[j], 3, els);

          for (size_t k=0; k < els.getSize(); ++k)
          {
            apf::MeshEntity* other_el = els[k];
            if (m_apf_data.m->toModel(other_el) != me)
            {
              ArrayNode node = getArrayNode(other_el, verts[j]);
              m_local_dest_nodes.push_back(node);
              end_idx++;
            }
          }

          if (end_idx != start_idx)
          {
            m_local_dest_node_ptrs.push_back(end_idx);
            m_local_src_nodes.push_back(getArrayNode(el, verts[j]));
          }
        } 
      }
    }
  }
}


DirichletUpdateMap::ArrayNode DirichletUpdateMap::getArrayNode(apf::MeshEntity* vert)
{
  apf::Adjacent els;
  m_apf_data.m->getAdjacent(vert, 3, els);
  apf::MeshEntity* el = els[0];

  return getArrayNode(el, vert);
}


void DirichletUpdateMap::getArrayNodes(apf::MeshEntity* vert, std::vector<ArrayNode>& nodes)
{
  apf::Adjacent els;
  m_apf_data.m->getAdjacent(vert, 3, els);

  for (auto el : els)
  {
    nodes.push_back(getArrayNode(el, vert));
  }
}

DirichletUpdateMap::ArrayNode DirichletUpdateMap::getArrayNode(apf::MeshEntity* el, apf::MeshEntity* vert)
{

  DofInt dof_num = apf::getNumber(m_apf_data.dof_nums, vert, 0, 0);
  
  int el_num = apf::getNumber(m_apf_data.el_nums, el, 0, 0);
  int el_num_local = m_mesh->m_elnums_global_to_local[el_num];
  int vol_group_num = apf::getNumber(m_apf_data.vol_groups, el, 0, 0);
  VolumeGroup& vol_group = m_mesh->getElements(vol_group_num);

  for (int i=0; i < vol_group.getNumSolPtsPerElement(); ++i)
    if (vol_group.nodenums[el_num_local][i] == dof_num)
      return ArrayNode{vol_group_num, el_num_local, i};

  throw std::runtime_error("unable to find array node for given vert"); 
}

bool DirichletUpdateMap::isSharedOrGhost(apf::MeshEntity* e)
{
  return m_apf_data.m->isShared(e) || m_apf_data.m->isGhosted(e) || m_apf_data.m->isGhost(e);
}


}