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

// Creates a ghost layer of elements.
// Specifically, add a layer of elements connect to local elements via a vertex.
// Some definitions:
//   An entity has isGhost() == true if it is created in this function
//   An entity has isGhosted() == true if it is created on another process in this
//   function.  For shared entitities that are ghosted, *all* the sharers have 
//   isGhosted() == true.
// The function getGhosts returns a map of (mpi_rank, apf::MeshEntity*) that contains
// the entities that are both isGhost() and isGhosted().  Therefore, this
// ghosting information can be used for symmetric communication patterns.
// getResidence() returns the set of procs the entity existed on *before*
// ghosting.  It does this for both regular and ghosted entities (ie.
// The ranks the entity was ghosted to are not included).
void createGhosting(apf::Mesh2* mesh, MPI_Comm comm);

}  // namespace

#endif