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
          if (mesh->isShared(e))
          {
            int owner_rank = mesh->getOwner(e);
            if (owner_rank == m_comm_rank)
            {
              std::cout << "entity " << e << " is shared and owned, sending to other ranks" << std::endl;
              copies.clear();
              mesh->getRemotes(e, copies);
              packCopies(e, copies);
              for (auto& p : copies)
                std::cout << "  sending to rank " << p.first << ", local entity " << p.second << std::endl;
            } else
            {
              std::cout << "entity " << e << " is shared and not owned, expecting to receive from rank " << m_field.getMesh()->getOwner(e) << std::endl;
              countReceives(e);
            }
          }

          if (mesh->isGhost(e))
          {
            std::cout << "entity " << e << " is ghost and not owned, expecting to receive from rank " << m_field.getMesh()->getOwner(e) << std::endl;

            countReceives(e);
          } else if (mesh->isGhosted(e) && mesh->getOwner(e) == m_comm_rank)  // not sure if the getOwner condition is required,
                                                                              // it might be implied by isGhosted
          {
            std::cout << "entity " << e << " is shared and owned, sending to other ranks" << std::endl;

            copies.clear();
            mesh->getGhosts(e, copies);
            packCopies(e, copies);

            for (auto& p : copies)
              std::cout << "  sending to rank " << p.first << ", local entity " << p.second << std::endl;
          }
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
            m_send_vals[rank].push_back(m_field(entity_local, i, c));
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
          int count = m_recv_entities[i].size() * sizeof(apf::MeshEntity*);
          MPI_Irecv(m_recv_entities[i].data(), count, MPI_BYTE, i, tag, m_comm, &(recv_reqs_entities[idx]));

          count = m_recv_vals[i].size() * sizeof(T);
          MPI_Irecv(m_recv_vals[i].data(), count, MPI_BYTE, i, tag, m_comm, &(recv_reqs_vals[idx++]));
        }

      idx = 0;
      for (int i=0; i < m_comm_size; ++i)
        if (m_send_entities[i].size() > 0)
        {
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
          auto& vals = m_recv_vals[i];
          int idx = 0;
          for (auto& e : m_recv_entities[i])
          {
            std::cout << "receiving entity " << e << std::endl;
            for (int i=0; i < countNodesOn(e); ++i)
              for (int c=0; c < m_field.getNumComponents(); ++c)
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