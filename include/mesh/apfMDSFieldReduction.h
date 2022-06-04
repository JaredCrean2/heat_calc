#ifndef APF_MDS_FIELD_REDUCTION_H
#define APF_MDS_FIELD_REDUCTION_H

#include "mesh/apfMDSField.h"
#include "mpi.h"

// for source file
#include "PCU.h"


namespace fast_field {

#ifdef MESH_USE_MDS_NUMBERING

template <typename T>
struct assign {
  constexpr T operator()(const T& lhs, const T& rhs) const { return lhs; }
};

template <typename T, typename MeshWrapper, typename Func=std::plus<T>>
class ApfMDSFieldReduction
{
  public:
    // apply func to parallel entities.  Can include ghost entities and/or shared entities
    explicit ApfMDSFieldReduction(ApfMDSField<T, MeshWrapper>& field, const Func& func=std::plus<T>(),
                                  bool include_shared=true, bool include_ghosts=true) :
      m_field(field),
      m_func(func),
      m_comm(PCU_Get_Comm()),
      m_include_shared(include_shared),
      m_include_ghosts(include_ghosts)
    {
      if (!include_shared && !include_ghosts)
        throw std::runtime_error("neither shared nor ghost entities selected for reduction.  Reduction would have no effect");
    
      MPI_Comm_rank(m_comm, &m_comm_rank);
      MPI_Comm_size(m_comm, &m_comm_size);
      m_send_entities.resize(m_comm_size);
      m_send_vals.resize(m_comm_size);
      m_recv_entity_counts.resize(m_comm_size, 0);
      m_recv_value_counts.resize(m_comm_size, 0);
      m_recv_entities.resize(m_comm_size);
      m_recv_vals.resize(m_comm_size);
    }

    void execute()
    {
      packBuffers();
      doCommunication();
    }

  private:

    void packBuffers()
    {
      apf::Mesh2* mesh = m_field.getMesh();
      apf::Copies copies, copies_ghost;

      for (int dim=0; dim < 3; ++dim)
      {
        if (!m_field.getFieldShape()->hasNodesIn(dim))
          continue;

        apf::MeshIterator* it = mesh->begin(dim);
        apf::MeshEntity* e;
        while ( (e = mesh->iterate(it)))
        {
          copies.clear();

          if (mesh->isShared(e) && m_include_shared)
            mesh->getRemotes(e, copies);

          if ( (mesh->isGhost(e) || mesh->isGhosted(e)) && m_include_ghosts)
          {
            copies_ghost.clear();
            mesh->getGhosts(e, copies_ghost);

            for (auto& p : copies_ghost)
              copies.insert(p);
          }
          
          packCopies(e, copies);
          countReceives(e, copies);
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

    void countReceives(apf::MeshEntity* entity_local, apf::Copies& copies)
    {
      int nvals = countNodesOn(entity_local) * m_field.getNumComponents();
      for (auto& p : copies)
      {
        int rank = p.first;
        m_recv_entity_counts[rank]++;
        m_recv_value_counts[rank] += nvals;
      }
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
            for (int i=0; i < countNodesOn(e); ++i)
              for (int c=0; c < m_field.getNumComponents(); ++c)
                m_field(e, i, c) = m_func(m_field(e, i, c), vals[idx++]);
        }
    }

    ApfMDSField<T, MeshWrapper>& m_field;
    const Func& m_func;
    MPI_Comm m_comm;
    bool m_include_shared;
    bool m_include_ghosts;
    int m_comm_rank;
    int m_comm_size;

    std::vector< std::vector<apf::MeshEntity*> > m_send_entities;
    std::vector< std::vector<T> > m_send_vals;
    std::vector<int> m_recv_entity_counts;
    std::vector<int> m_recv_value_counts;
    std::vector< std::vector<apf::MeshEntity*> > m_recv_entities;
    std::vector<std::vector<T> > m_recv_vals;
};

template <typename T, typename MeshWrapper, typename Func>
ApfMDSFieldReduction<T, MeshWrapper, Func> make_field_reduction(ApfMDSField<T, MeshWrapper>& field, const Func& func)
{
  return ApfMDSFieldReduction<T, MeshWrapper, Func>(field, func);
}

} // namespace

#endif

#endif