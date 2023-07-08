#ifndef HEAT_CALC_UTILS_PARALLEL_EXCHANGE_H
#define HEAT_CALC_UTILS_PARALLEL_EXCHANGE_H

#include "mpi.h"
#include "mpi_utils.h"
#include <iostream>
#include <vector>

template <typename T>
class ParallelExchange
{
  public:
    ParallelExchange(MPI_Comm comm, int tag) :
      m_comm(comm),
      m_tag(tag),
      m_send_bufs(commSize(comm)),
      m_recv_bufs(commSize(comm)),
      m_send_reqs(commSize(comm), MPI_REQUEST_NULL),
      m_recv_reqs(commSize(comm), MPI_REQUEST_NULL)
    {}

    std::vector<T>& getSendBuf(int rank) { return m_send_bufs[rank]; }

    const std::vector<T>& getSendBuf(int rank) const { return m_send_bufs[rank]; }

    std::vector<T>& getRecvBuf(int rank) { return m_recv_bufs[rank]; }

    const std::vector<T>& getRecvBuf(int rank) const { return m_recv_bufs[rank]; }

    void startCommunication()
    {
      for (size_t i=0; i < m_recv_bufs.size(); ++i)
        if (m_recv_bufs[i].size() > 0)
          MPI_Irecv(m_recv_bufs[i].data(), m_recv_bufs[i].size() * sizeof(T), MPI_BYTE, i, m_tag, m_comm, &(m_recv_reqs[i]));

      for (size_t i=0; i < m_send_bufs.size(); ++i)
        if (m_send_bufs[i].size() > 0)
          MPI_Isend(m_send_bufs[i].data(), m_send_bufs[i].size() * sizeof(T), MPI_BYTE, i, m_tag, m_comm, &(m_send_reqs[i]));
    }

    template <typename Func>
    void finishCommunication(Func func)
    {
      int numRecvs = 0;
      for (size_t i=0; i < m_recv_bufs.size(); ++i)
        if (m_recv_bufs[i].size() > 0)
          numRecvs++;

      std::cout << "numRecvs = " << numRecvs << std::endl;

      for (int i=0; i < numRecvs; ++i)
      {
        std::cout << "i = " << i << std::endl;
        int rank;
        MPI_Waitany(m_recv_reqs.size(), m_recv_reqs.data(), &rank, MPI_STATUS_IGNORE);
        std::cout << "rank = " << rank << std::endl;
        func(rank, m_recv_bufs[rank]);
      }

      MPI_Waitall(m_send_reqs.size(), m_send_reqs.data(), MPI_STATUSES_IGNORE);
    }

  private:
    MPI_Comm m_comm;
    int m_tag;

    std::vector<std::vector<T>> m_send_bufs;
    std::vector<std::vector<T>> m_recv_bufs;
    std::vector<MPI_Request> m_send_reqs;
    std::vector<MPI_Request> m_recv_reqs;

};

#endif