#ifndef APF_MDS_FIELD_SYNCHRONIZE_H
#define APF_MDS_FIELD_SYNCHRONIZE_H

#include "mesh/apfMDSField.h"
#include "mpi.h"

// for source file
#include "PCU.h"


namespace fast_field {

#ifdef MESH_USE_MDS_NUMBERING


template <typename T, typename MeshWrapper>
class ApfMDSFieldSyncrhonizer
{
  public:
    explicit ApfMDSFieldSyncrhonizer(ApfMDSField<T, MeshWrapper>& field) :
      m_field(field),
      m_comm(PCU_Get_Comm())
    {
      MPI_Comm_rank(m_comm, &m_comm_rank);
      MPI_Comm_size(m_comm, &m_comm_size);
      m_send_entities.resize(m_comm_size);
      m_send_vals.resize(m_comm_size);
      m_recv_entity_counts.resize(m_comm_size, 0);
      m_recv_value_counts.resize(m_comm_size, 0);
      m_recv_entities.resize(m_comm_size);
      m_recv_vals.resize(m_comm_size);
    }

    void synchronize()
    {
      std::cout << "\nsynchronizing field" << m_field.getName() <<  std::endl;
      packBuffers();
      doCommunication();
    }

  private:

    void packBuffers()
    {
      apf::Mesh2* mesh = m_field.getMesh();
      apf::Copies copies;

      for (int dim=0; dim <= 3; ++dim)
      {
        if (!m_field.getFieldShape()->hasNodesIn(dim))
          continue;

        apf::MeshIterator* it = mesh->begin(dim);
        apf::MeshEntity* e;
        while ( (e = mesh->iterate(it)))
        {
          bool isParallel = mesh->isShared(e) || mesh->isGhost(e) || mesh->isGhosted(e);
          if (isParallel)
          {
            bool isOwned = mesh->getOwner(e) == m_comm_rank;
            if (isOwned)
            {
              copies.clear();
              mesh->getRemotes(e, copies);
              mesh->getGhosts(e, copies);
              packCopies(e, copies);
            } else
            {
              countReceives(e);
            }


          }
/*
          if (mesh->isShared(e))
          {
            std::cout << "packing shares for entity " << e << std::endl;
            int owner_rank = mesh->getOwner(e);
            if (owner_rank == m_comm_rank)
            {
              copies.clear();
              mesh->getRemotes(e, copies);
              packCopies(e, copies);
            } else
            {
              countReceives(e);
            }
          }

          // This problem is this: our definitions are that
          //  1) getGhosts on the owner returns the ghosts, and getGhosts on the ghosts returns (at least) the owner
          //  2) getGhosts should return the same values on the sharers as on the ghosts (ie. everyone knows about everyone)
          // For more than 2 procs, this implies that two shared entities can also be ghost entities if the entity is ghosted
          // to a third process.  It happens during the updateGhosts() function: the ghost initially has the owner is its
          // getGhosts list, it then sends it getGhosts list to the other ghosts, which includes the shared entities.
          // So now the shared entities has the owner in its getGhosts list.  One more round results in the owner having
          // The shared entity in the ghost list.
          if (mesh->isGhost(e))  
          {
            std::cout << "counting ghost receives for entity " << e << std::endl;
            countReceives(e);
          } else if (mesh->isGhosted(e) && mesh->getOwner(e) == m_comm_rank)  // not sure if the getOwner condition is required,
                                                                              // it might be implied by isGhosted
          {
            std::cout << "packing ghosts for entity " << e << std::endl;
            copies.clear();
            mesh->getGhosts(e, copies);
            packCopies(e, copies);
          }
*/          
        }
        mesh->end(it);
      }
    }

    void packCopies(apf::MeshEntity* entity_local, apf::Copies& copies)
    {

      for (auto& p : copies)
      {
        int rank = p.first;

        m_send_entities[rank].push_back(p.second);
        for (int i=0; i < countNodesOn(entity_local); ++i)
          for (int c=0; c < m_field.getNumComponents(); ++c)
          {
            //if (entity_local == reinterpret_cast<apf::MeshEntity*>((void*)0x1f9))
            //  std::cout << "sending value " << m_field(entity_local, i, c) << " to rank " << rank << " for entity " << p.second << " from local entity " << entity_local << std::endl;
            m_send_vals[rank].push_back(m_field(entity_local, i, c));
          }
      }
    }

    void countReceives(apf::MeshEntity* entity_local)
    {
      int owner_rank = m_field.getMesh()->getOwner(entity_local);
      m_recv_entity_counts[owner_rank]++;
      m_recv_value_counts[owner_rank] += countNodesOn(entity_local) * m_field.getNumComponents();
    }

    int countNodesOn(apf::MeshEntity* e)
    {
      auto etype = m_field.getMesh()->getType(e);
      return m_field.getFieldShape()->countNodesOn(etype);
    }

    void doCommunication()
    {
      int num_sends = countSends();
      int num_recvs = sizeReceiveBuffers();
      std::vector<MPI_Request> send_reqs_entities(num_sends), send_reqs_vals(num_sends);
      std::vector<MPI_Request> recv_reqs_entities(num_recvs), recv_reqs_vals(num_recvs);
      const int tag = 10100;


      int idx = 0;
      for (int i=0; i < m_comm_size; ++i)
        if (m_recv_entity_counts[i] > 0)
        {
          std::cout << "expecting to receive " << m_recv_entities[i].size() << " entities and " << m_recv_vals[i].size() << " values from rank " << i << std::endl;
          int count = m_recv_entities[i].size() * sizeof(apf::MeshEntity*);
          MPI_Irecv(m_recv_entities[i].data(), count, MPI_BYTE, i, tag, m_comm, &(recv_reqs_entities[idx]));

          count = m_recv_vals[i].size() * sizeof(T);
          MPI_Irecv(m_recv_vals[i].data(), count, MPI_BYTE, i, tag, m_comm, &(recv_reqs_vals[idx++]));
        }

      idx = 0;
      for (int i=0; i < m_comm_size; ++i)
        if (m_send_entities[i].size() > 0)
        {
          std::cout << "sending " << m_send_entities[i].size() << " entities and " << m_send_vals[i].size() << " values from rank " << i << std::endl;

          int count = m_send_entities[i].size() * sizeof(apf::MeshEntity*);
          MPI_Isend(m_send_entities[i].data(), count, MPI_BYTE, i, tag, m_comm, &(send_reqs_entities[idx]));

          count = m_send_vals[i].size() * sizeof(T);
          MPI_Isend(m_send_vals[i].data(), count, MPI_BYTE, i, tag, m_comm, &(send_reqs_vals[idx++]));
        }

      MPI_Waitall(num_recvs, recv_reqs_entities.data(), MPI_STATUSES_IGNORE);
      MPI_Waitall(num_recvs, recv_reqs_vals.data(),     MPI_STATUSES_IGNORE);
      unpackRecvBuffers();
      MPI_Waitall(num_sends, send_reqs_entities.data(), MPI_STATUSES_IGNORE);
      MPI_Waitall(num_sends, send_reqs_vals.data(),     MPI_STATUSES_IGNORE);
    }

    int countSends()
    {
      auto f = [](const std::vector<apf::MeshEntity*>& v) { return v.size() > 0; };
      return std::count_if(m_send_entities.begin(), m_send_entities.end(), f);
    }

    int sizeReceiveBuffers()
    {
      int num_recvs = 0;
      for (int i=0; i < m_comm_size; ++i)
        if (m_recv_entity_counts[i] > 0)
        {
          m_recv_entities[i].resize(m_recv_entity_counts[i]);
          m_recv_vals[i].resize(m_recv_value_counts[i]);
          num_recvs++;
        }
        
      return num_recvs;
    }

    void unpackRecvBuffers()
    {
      for (int i=0; i < m_comm_size; ++i)
        if (m_recv_entities[i].size() > 0)
        {
          std::cout << "\nreceiving from rank " << i << std::endl;
          auto& vals = m_recv_vals[i];
          int idx = 0;
          for (auto& e : m_recv_entities[i])
            for (int i=0; i < countNodesOn(e); ++i)
              for (int c=0; c < m_field.getNumComponents(); ++c)
              {
                std::cout << "receiving val " << vals[idx] << " for entity " << e << std::endl;
                m_field(e, i, c) = vals[idx++];
              }
        }
    }

    ApfMDSField<T, MeshWrapper>& m_field;
    MPI_Comm m_comm;
    int m_comm_rank;
    int m_comm_size;
    std::vector< std::vector<apf::MeshEntity*> > m_send_entities;
    std::vector< std::vector<T> > m_send_vals;
    std::vector<int> m_recv_entity_counts;
    std::vector<int> m_recv_value_counts;
    std::vector< std::vector<apf::MeshEntity*> > m_recv_entities;
    std::vector<std::vector<T> > m_recv_vals;
};

template <typename T, typename MeshWrapper>
ApfMDSFieldSyncrhonizer<T, MeshWrapper> make_field_synchronizer(ApfMDSField<T, MeshWrapper>& field)
{
  return ApfMDSFieldSyncrhonizer<T, MeshWrapper>(field);
}

} // namespace

#endif

#endif