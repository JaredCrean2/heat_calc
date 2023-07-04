#ifndef HEAT_CALC_MESH_GHOSTING_H
#define HEAT_CALC_MESH_GHOSTING_H

#include "apfMesh2.h"
#include "error_handling.h"
#include "mpi.h"
#include "PCU.h"

// Usage of PCU:
// PCU_Comm_Begin()  // Start communication phase
// PCU_Comm_Pack(dest, data)  // pack data into buffers
// PCU_Comm_Send();  // send the data
// while (PCU_Comm_Listen())
// {
//   // now PCU has received a buffer and set it to be the active buffer
//   PCU_Comm_Sender()  // rank of sender
//   while (!PCU_Buffer_Unpacked())
//   {
//     PCU_Comm_Unpack(void* data, size_t size) // unpacks data from the active buffer
//   }                                          // and puts it into data
//
//  }

template <typename T>
void PCU_Comm_Pack(int dest_rank, const T& data)
{
  int status = PCU_Comm_Pack(dest_rank, reinterpret_cast<const void*>(&data), sizeof(T));
  assertAlways(status == PCU_SUCCESS, "PCU_Comm_Pack failed");
}

template <typename T>
T PCU_Comm_Unpack()
{
  T val;
  PCU_Comm_Unpack(static_cast<void*>(&val), sizeof(T));
  return val;
}

namespace Mesh
{

void createGhosting(apf::Mesh2* mesh, MPI_Comm comm);

}  // namespace

#endif