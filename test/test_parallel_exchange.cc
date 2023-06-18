#include "gtest/gtest.h"
#include "utils/parallel_exchange.h"

namespace {

int get_value(int send_rank, int dest_rank, int idx)
{
  return send_rank + 2*dest_rank + 3*idx;
}

}

TEST(ParallelExchange, AlltoAll)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  ParallelExchange<int> exchanger(comm, 1000);

  int myrank = commRank(comm);
  for (int dest_rank=0; dest_rank < commSize(comm); ++dest_rank)
  {
    int buf_size = myrank + dest_rank;
    for (int i=0; i < buf_size; ++i)
      exchanger.getSendBuf(dest_rank).push_back(get_value(myrank, dest_rank, i));

  }

  for (int send_rank=0; send_rank < commSize(comm); ++send_rank)
    exchanger.getRecvBuf(send_rank).resize(send_rank + myrank);

  exchanger.startCommunication();

  auto f = [&](int rank, const std::vector<int>& buf)
  {
    EXPECT_EQ(buf.size(), size_t(rank + myrank));

    for (int i=0; i < (rank + myrank); ++i)
      EXPECT_EQ(buf[i], get_value(rank, myrank, i));
  };

  exchanger.finishCommunication(f);

}