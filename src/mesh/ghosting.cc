#include "apf.h"
#include "apfMesh2.h"
#include "error_handling.h"
#include <PCU.h>
#include <algorithm>
#include <apfMesh.h>
#include <cassert>
#include "mpi_utils.h"
#include <array>
#include <iostream>

#include "mesh/ghosting.h"


namespace Mesh {

apf::Vector3 computeCentroid(apf::Mesh2* mesh, apf::MeshEntity* e)
{
  if (mesh->getType(e) == apf::Mesh::Type::VERTEX)
  {
    apf::Vector3 pt;
    mesh->getPoint(e, 0, pt);
    return pt;
  } else
  {
    apf::Vector3 centroid;
    apf::Downward down;
    int ndown = mesh->getDownward(e, 0, down);
    for (int i=0; i < ndown; ++i)
    {
      centroid += computeCentroid(mesh, down[i]);
    }

    return centroid / ndown;
  }
}

class GhostingCreator
{
  public:
    GhostingCreator(apf::Mesh2* mesh, MPI_Comm comm) :
      m_mesh(mesh),
      m_comm(comm)
    {
      m_is_ghost_tag = m_mesh->findTag("ghost_tag");
      m_is_ghosted_tag = m_mesh->findTag("ghosted_tag");

      if (!m_is_ghost_tag)
        m_is_ghost_tag = m_mesh->createIntTag("ghost_tag",1);

      if (!m_is_ghosted_tag)
        m_is_ghosted_tag = m_mesh->createIntTag("ghosted_tag",1);      
    }

    void create()
    {
      std::vector<apf::MeshEntity*> els_to_ghost;
      getGhostElements(els_to_ghost);
      for (auto& el : els_to_ghost)
      {
        std::cout << "ghosting el " << el << " with centroid " << computeCentroid(m_mesh, el) << std::endl;
      }

      //std::cout << "\nGhosting verts" << std::endl;
      //ghostVerts(els_to_ghost);  

      std::cout << "\nghosting higher dimension entities" << std::endl;
      ghostEntities(els_to_ghost);
      m_mesh->acceptChanges();
    }

    void getGhostElements(std::vector<apf::MeshEntity*>& els)
    {
      // loop over verts
      //   if isShared
      //     get all upward elements
      //     push_back into vector
      apf::MeshIterator* mit = m_mesh->begin(0);
      apf::MeshEntity* vert;
      apf::Adjacent vert_els;

      while ( (vert = m_mesh->iterate(mit)) )
      {
        if (m_mesh->isShared(vert))
        {
          vert_els.resize(0);
          m_mesh->getAdjacent(vert, 3, vert_els);
          for (int i=0; i < vert_els.getSize(); ++i)
            els.push_back(vert_els[i]);
        }
      }
      m_mesh->end(mit);

      std::sort(els.begin(), els.end());
      auto it = std::unique(els.begin(), els.end());
      els.erase(it, els.end());
    }


    void getEntityDestinations(const std::vector<apf::MeshEntity*>& els_to_ghost, int dim,
                             std::map<apf::MeshEntity*, std::vector<int>>& entity_dests)
    {
      apf::Copies shared_entities;
      for (auto& el : els_to_ghost)
      {
        apf::Downward verts, entities;
        int nverts = m_mesh->getDownward(el, 0, verts);
        int nentities = get_downward_or_self(el, dim, entities);

        for (int i=0; i < nverts; ++i)
          if (m_mesh->isShared(verts[i]))
          {
            shared_entities.clear();
            m_mesh->getRemotes(verts[i], shared_entities);
            for (int j=0; j < nentities; ++j)
              for (auto& p : shared_entities)
              {
                entity_dests[entities[j]].push_back(p.first);
              }
          }
      }

      sort_and_unique(entity_dests);
    }

    int get_downward_or_self(apf::MeshEntity* e, int dim, apf::MeshEntity** down)
    {
      if (apf::Mesh::typeDimension[m_mesh->getType(e)] == dim)
      {
        down[0] = e;
        return 1;
      } else
      {
        return m_mesh->getDownward(e, dim, down);
      }
    }

    void sort_and_unique(std::map<apf::MeshEntity*, std::vector<int>>& entity_dests)
    {
      for (auto & p : entity_dests)
      {
        std::vector<int>& dests = p.second;
        std::sort(dests.begin(), dests.end());
        auto it = std::unique(dests.begin(), dests.end());
        dests.erase(it, dests.end());
      }
    }

    void getEntityDestinationsOnOwnerRank(std::map<apf::MeshEntity*, std::vector<int>>& entity_dests)
    {
      apf::Copies shared_entities;
      PCU_Comm_Begin();
      for (auto it = entity_dests.begin(); it != entity_dests.end(); )
      {
        apf::MeshEntity* entity      = it->first;
        std::vector<int>& dest_ranks = it->second;

        if (!m_mesh->isOwned(entity))
        {
          int owner_rank = m_mesh->getOwner(entity);

          shared_entities.clear();
          m_mesh->getRemotes(entity, shared_entities);
          apf::MeshEntity* entity_on_owner = shared_entities.at(owner_rank);

          PCU_Comm_Pack(owner_rank, entity_on_owner);
          int ndests = dest_ranks.size();
          PCU_Comm_Pack(owner_rank, ndests);
          for (auto& rank : dest_ranks)
            PCU_Comm_Pack(owner_rank, rank);

          it = entity_dests.erase(it);
        } else
        {
          it++;
        }
      }

      PCU_Comm_Send();

      int myrank = commRank(m_comm);
      while (PCU_Comm_Listen())
      {
        while (!PCU_Comm_Unpacked())
        {
          auto entity = PCU_Comm_Unpack<apf::MeshEntity*>();
          std::vector<int>& dest_ranks = entity_dests[entity];

          auto ndests = PCU_Comm_Unpack<int>();
          for (int i=0; i < ndests; ++i)
          {
            int dest_rank = PCU_Comm_Unpack<int>();
            if (dest_rank != myrank)
              dest_ranks.push_back(dest_rank);
          }
        }
      }

      sort_and_unique(entity_dests);
    }


    void packVertex(apf::MeshEntity* vert, int dest, const apf::Copies& shared_entities)
    {
      // pack coordinates
      // if already exists, pack remote, else nullptr  //No: if the vert already exists, it is a share, not a ghost
      // pack vert
      // if vert is shared
      //   pack remotes
      // set isGhosted flag on vert

      apf::Vector3 pt;
      m_mesh->getPoint(vert, 0, pt);
      std::array<double, 3> pt_array = {pt.x(), pt.y(), pt.z()};
      apf::ModelEntity* me = m_mesh->toModel(vert);
      int me_dim = m_mesh->getModelType(me);
      int me_tag = m_mesh->getModelTag(me);

      m_mesh->getPoint(vert, 0, pt);

      PCU_Comm_Pack(dest, pt_array);
      PCU_Comm_Pack(dest, me_dim);
      PCU_Comm_Pack(dest, me_tag);
      PCU_Comm_Pack(dest, vert);
      int nshared = shared_entities.size();
      PCU_Comm_Pack(dest, nshared);
      for (auto& p : shared_entities)
      {
        PCU_Comm_Pack(dest, p.first);
        PCU_Comm_Pack(dest, p.second);
      }

      // need to set isGhosted flag on all shared, but do that on the reverse
      set_is_ghosted(vert);
    }

    apf::MeshEntity* unpackVertex()
    {
      // unpack coords
      // unpack existing vert
      // create vertex if no existing vert
      // set isGhost flag
      // addGhost(sendingRank, sent vert)
      // unpack remotes, addGhost

      auto pt_array = PCU_Comm_Unpack<std::array<double, 3>>();
      int me_dim = PCU_Comm_Unpack<int>();
      int me_tag = PCU_Comm_Unpack<int>();
      int sender_rank = PCU_Comm_Sender();

      apf::MeshEntity* vert_sender = PCU_Comm_Unpack<apf::MeshEntity*>();
      apf::Vector3 pt(pt_array[0], pt_array[1], pt_array[2]);
      apf::ModelEntity* me = m_mesh->findModelEntity(me_dim, me_tag);

      apf::MeshEntity* vert = m_mesh->createVert(me);
      m_mesh->setPoint(vert, 0, pt);
      set_is_ghost(vert);
      std::cout << "creating new vert " << vert << " at " << pt << std::endl;


      m_mesh->addGhost(vert, sender_rank, vert_sender);
      int nshared = PCU_Comm_Unpack<int>();
      for (int i=0; i < nshared; ++i)
      {
        int rank = PCU_Comm_Unpack<int>();
        apf::MeshEntity* vert_shared = PCU_Comm_Unpack<apf::MeshEntity*>();
        m_mesh->addGhost(vert, rank, vert_shared);
      }

      return vert;
      // pack newly-created vert back to sender
    }

    void returnGhostsToOwner(const std::vector<std::vector<apf::MeshEntity*>>& sent_entities,
                             const std::vector<std::vector<apf::MeshEntity*>>& recv_entities)
    {
      PCU_Comm_Begin();
      for (int rank=0; rank < int(recv_entities.size()); ++rank)
        for (auto& entity : recv_entities[rank])
        {
          PCU_Comm_Pack(rank, entity);
        }

      PCU_Comm_Send();

      while (PCU_Comm_Listen())
      {
        int sender_rank = PCU_Comm_Sender();
        int idx = 0;
        while (!PCU_Comm_Unpacked())
        {
          auto entity_local  = sent_entities[sender_rank][idx++];
          auto entity_remote = PCU_Comm_Unpack<apf::MeshEntity*>();
          m_mesh->addGhost(entity_local, sender_rank, entity_remote);
        }
      }
    }


    int updateGhosts(int dim)
    {
      // all ghosts send  getGhosts to all ghosts (plus self)
      // all ghosts receive and add new ghosts to getGhosts
      // if !isGhost, set isGhosted
      // return number of new ghosts
      // call this function as many times as necessary until number of new ghosts is zero

      int myrank = commRank(m_comm);
      PCU_Comm_Begin();
      apf::Copies ghosts;
      apf::MeshIterator* it = m_mesh->begin(dim);
      apf::MeshEntity* entity;
      while ( (entity = m_mesh->iterate(it)) )
      {
        if (m_mesh->isGhost(entity) || m_mesh->isGhosted(entity))
        {
          ghosts.clear();
          m_mesh->getGhosts(entity, ghosts);
          ghosts[myrank] = entity;

          for  (auto& p : ghosts)
          {
            int dest_rank = p.first;
            if (dest_rank != myrank)
            {
              for (auto& p2 : ghosts)
              {
                apf::MeshEntity* entity_remote = p.second;
                int rank_ghost = p2.first;
                apf::MeshEntity* entity_ghost = p2.second;
                PCU_Comm_Pack(dest_rank, entity_remote);
                PCU_Comm_Pack(dest_rank, rank_ghost);
                PCU_Comm_Pack(dest_rank, entity_ghost);
              }
            }
          }
        }
      }
      m_mesh->end(it);

      PCU_Comm_Send();

      int num_new_ghosts = 0;
      while (PCU_Comm_Listen())
      {
        while (!PCU_Comm_Unpacked())
        {
          auto entity_local = PCU_Comm_Unpack<apf::MeshEntity*>();
          int rank_ghost    = PCU_Comm_Unpack<int>();
          auto entity_ghost = PCU_Comm_Unpack<apf::MeshEntity*>();

          ghosts.clear();
          m_mesh->getGhosts(entity_local, ghosts);
          if (ghosts.count(rank_ghost) == 0)
          {
            m_mesh->addGhost(entity_local, rank_ghost, entity_ghost);
            num_new_ghosts++;
            if (!m_mesh->isGhost(entity_local) && !m_mesh->isGhosted(entity_local))
            {
              std::cout << "setting entity " << entity_local << " to is_ghosted from updateGhosts" << std::endl;
              set_is_ghosted(entity_local);
            }
          }
        }
      }

      int num_new_ghosts_global = 0;
      MPI_Allreduce(&num_new_ghosts, &num_new_ghosts_global, 1, MPI_INT, MPI_SUM, m_comm);

      return num_new_ghosts_global;
    }

    void ghostEntities(const std::vector<apf::MeshEntity*>& els_to_ghost)
    {
      for (int dim=0; dim <= 3; ++dim)
      {
        std::map<apf::MeshEntity*, std::vector<int>> entity_dests;
        getEntityDestinations(els_to_ghost, dim, entity_dests);
        getEntityDestinationsOnOwnerRank(entity_dests);

        createEntitiesFromOwner(entity_dests, dim);

        int num_new_ghosts = 0;
        do
        {
          num_new_ghosts = updateGhosts(dim);
        } while (num_new_ghosts > 0);
      }

    }

    void createEntitiesFromOwner(std::map<apf::MeshEntity*, std::vector<int>>& entity_dests, int dim)
    {
      // get vector of edges
      // get destinations for each edge from getGhosts(verts)
      // send all destinations to owner of edge
      // (on owner) for each edge 
      //   pack edge to destination
      // set isGhosted

      PCU_Comm_Begin();
      std::vector<std::vector<apf::MeshEntity*>> sent_entities(commSize(m_comm));
      apf::Copies shared_entities;
      for (auto& p : entity_dests)
      {
        apf::MeshEntity* entity      = p.first;
        std::vector<int>& dest_ranks = p.second;

        shared_entities.clear();
        m_mesh->getRemotes(entity, shared_entities);
        for (int dest_rank : dest_ranks)
        {
          if (shared_entities.count(dest_rank) == 0)
          {
            packEntity(entity, dest_rank, shared_entities);
            sent_entities[dest_rank].push_back(entity);
          }
        }
      }

      PCU_Comm_Send();

      std::vector<std::vector<apf::MeshEntity*>> received_entities(commSize(m_comm));
      while (PCU_Comm_Listen())
      {
        int sender_rank = PCU_Comm_Sender();
        while (!PCU_Comm_Unpacked())
        {
          received_entities[sender_rank].push_back(unpackEntity(dim));
        }
      }

      returnGhostsToOwner(sent_entities, received_entities);
    }

    void packNonVertexEntity(apf::MeshEntity* entity, int dest, apf::Copies& entity_remotes)
    {
      std::cout << "sending entity " << entity << " at " << computeCentroid(m_mesh, entity) << " to rank " << dest << std::endl;
      // pack bounding verts on dest rank (from getRemotes/getGhosts of verts)
      // pack edge on remote (if present), or nullptr
      // pack remotes of edge

      int entity_type = m_mesh->getType(entity);
      int dim = apf::Mesh::typeDimension[entity_type];

      apf::Downward down_entities, down_entities_remote;
      int ndown = m_mesh->getDownward(entity, dim-1, down_entities);

      apf::Copies down_entity_remotes;
      for (int i=0; i < ndown; ++i)
      {
        std::cout << "down entity = " << down_entities[i] << " at " << computeCentroid(m_mesh, down_entities[i]) << std::endl;
        down_entity_remotes.clear();
        m_mesh->getRemotes(down_entities[i], down_entity_remotes);
        m_mesh->getGhosts(down_entities[i], down_entity_remotes);

        for (auto& p : down_entity_remotes)
          std::cout << "  rank " << p.first << ", entity " << p.second << std::endl;

        down_entities_remote[i] = down_entity_remotes.at(dest);
      }


      //apf::MeshEntity* vert1_dest = vert1_remotes.at(dest);
      //apf::MeshEntity* vert2_dest = vert2_remotes.at(dest);

      apf::ModelEntity* me = m_mesh->toModel(entity);
      int me_dim = m_mesh->getModelType(me);
      int me_tag = m_mesh->getModelTag(me);

      //apf::Copies entity_remotes;
      //m_mesh->getRemotes(entity, entity_remotes);

     
      PCU_Comm_Pack(dest, entity_type);
      for (int i=0; i < ndown; ++i)
        PCU_Comm_Pack(dest, down_entities_remote[i]);
        
      //PCU_Comm_Pack(dest, vert1_dest);
      //PCU_Comm_Pack(dest, vert2_dest);
      PCU_Comm_Pack(dest, me_dim);
      PCU_Comm_Pack(dest, me_tag);

      PCU_Comm_Pack(dest, entity);
      int nshared = entity_remotes.size();
      PCU_Comm_Pack(dest, nshared);
      for (auto& p : entity_remotes)
      {
        PCU_Comm_Pack(dest, p.first);
        PCU_Comm_Pack(dest, p.second);
      }

      set_is_ghosted(entity);
    }

    apf::MeshEntity* unpackNonVertexEntity()
    {
      // unpack 2 verts
      // unpack local edge ptr
      // if nullptr,
      //   create new edge
      //   set isGhosted
      // add all remotes as ghosts
      apf::Downward down_entities;

      int entity_type = PCU_Comm_Unpack<int>();
      int dim = apf::Mesh::typeDimension[entity_type];
      int ndown = apf::Mesh::adjacentCount[entity_type][dim-1];
      for (int i=0; i < ndown; ++i)
        down_entities[i] = PCU_Comm_Unpack<apf::MeshEntity*>();

      //auto vert1       = PCU_Comm_Unpack<apf::MeshEntity*>();
      //auto vert2       = PCU_Comm_Unpack<apf::MeshEntity*>();
      int me_dim       = PCU_Comm_Unpack<int>();
      int me_tag       = PCU_Comm_Unpack<int>();
      auto entity_remote = PCU_Comm_Unpack<apf::MeshEntity*>();
      int sender_rank  = PCU_Comm_Sender();

      apf::ModelEntity* me = m_mesh->findModelEntity(me_dim, me_tag);
      apf::MeshEntity* entity = m_mesh->createEntity(entity_type, me, down_entities);
      set_is_ghost(entity);

      std::cout << "created entity " << entity << " from rank " << sender_rank << ", entity " << entity_remote << std::endl;

      m_mesh->addGhost(entity, sender_rank, entity_remote);
      int nshared = PCU_Comm_Unpack<int>();
      for (int i=0; i < nshared; ++i)
      {
        int rank = PCU_Comm_Unpack<int>();
        auto remote_edge = PCU_Comm_Unpack<apf::MeshEntity*>();
        std::cout << " adding ghost " << rank << ", " << remote_edge << std::endl;
        m_mesh->addGhost(entity, rank, remote_edge);
      }

      return entity;
    }

    void packEntity(apf::MeshEntity* entity, int dest, apf::Copies& entity_remotes)
    {
      if (m_mesh->getType(entity) == apf::Mesh::Type::VERTEX)
        packVertex(entity, dest, entity_remotes);
      else
         packNonVertexEntity(entity, dest, entity_remotes);
    }
    
    apf::MeshEntity* unpackEntity(int dim)
    {
      if (dim == 0)
        return unpackVertex();
      else
        return unpackNonVertexEntity();
    }

    void set_is_ghost(apf::MeshEntity* vert)
    {
      int val = true;  // the value doesn't really matter, all that matters is that hasTag() returns true
      m_mesh->setIntTag(vert, m_is_ghost_tag, &val);
    }
    
    void set_is_ghosted(apf::MeshEntity* vert)
    {
      int val = true;  // the value doesn't really matter, all that matters is that hasTag() returns true
      std::cout << "setting entity " << vert << " to is_ghosted" << std::endl;
      m_mesh->setIntTag(vert, m_is_ghosted_tag, &val);
    }

  private:
    apf::Mesh2* m_mesh;
    MPI_Comm m_comm;

    apf::MeshTag* m_is_ghost_tag;
    apf::MeshTag* m_is_ghosted_tag;
};


void createGhosting(apf::Mesh2* mesh, MPI_Comm comm)
{
  GhostingCreator creator(mesh, comm);
  creator.create();
}

}