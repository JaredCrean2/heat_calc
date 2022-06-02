#ifndef MESH_APF_SPLIT
#define MESH_APF_SPLIT

#include "apfMesh2.h"
#include "error_handling.h"
#include "mpi.h"
#include <PCU.h>
#include <apfMDS.h>
#include <memory>
#include "apfZoltan.h"
#include "parma.h"

// this is an in-memory version of the zsplit.cc test

namespace mesh
{

//TODO: move function defintions to file
class ApfSplitter
{
  public:
    // mesh should be the mesh on the procs that the mesh is currently defined on,
    // and nullptr on the others
    // model should be the GMI model, valid on *all* procs
    explicit ApfSplitter(apf::Mesh2* mesh, gmi_model* model) :
      m_mesh(mesh),
      m_model(model)
    {}

    // on exit, the mesh will be split onto to_procs and PCU_Comm will be
    // a communicator with size to_procs
    // Users can optionally supply a communicator to use, otherwise a new
    // one will be created.
    // The returned mesh is the split mesh, valid on the new comm
    apf::Mesh2* split(int to_procs, MPI_Comm comm_new = MPI_COMM_NULL)
    {
      setupComm(comm_new, to_procs);
      int partition_factor = to_procs / PCU_Comm_Peers();

      apf::Migration* plan = nullptr;
      if (m_mesh)
        plan = getPlan(to_procs);

      PCU_Switch_Comm(m_comm_new);
      PCU_Barrier();
      // it is very unclear why repeatMdsMesh is called rather than expandMdsMesh and
      // then Mesh::migrate
      m_mesh = apf::repeatMdsMesh(m_mesh, m_model, plan, partition_factor);
      // repeatMdsMesh internally deletes the plan

      return m_mesh;
    }




  private:

    apf::Migration* getPlan(int to_procs)
    {
      int partition_factor = to_procs / PCU_Comm_Peers();
      apf::MeshTag* weights = Parma_WeighByMemory(m_mesh);
      std::shared_ptr<apf::Splitter> splitter(apf::makeZoltanSplitter(m_mesh, apf::GRAPH, apf::PARTITION, false));
      apf::Migration* plan = splitter->split(weights, 1.05, partition_factor);

      apf::removeTagFromDimension(m_mesh, weights, m_mesh->getDimension());
      m_mesh->destroyTag(weights);

      return plan;
    }

    void setupComm(MPI_Comm comm_new, int to_procs)
    {
      MPI_Comm_size(MPI_COMM_WORLD, &m_comm_size_world);
      MPI_Comm_rank(MPI_COMM_WORLD, &m_myrank_world);

      assertAlways(to_procs <= m_comm_size_world, "Cannot split mesh to more procs than size of MPI_COMM_WORLD");
      assertAlways(to_procs % PCU_Comm_Peers() == 0, "Number of procs to split mesh into must be a multiple of current size");

      if (comm_new == MPI_COMM_NULL)
      {    
        //TODO: in the multi-node case it would be better to distribute the
        //      mesh among nodes, but that would require additional calculation
        //      for the case where the new size is not divizible by the the comm size
        m_is_new = m_myrank_world < to_procs;
        if (to_procs == m_comm_size_world)
          m_comm_new = MPI_COMM_WORLD;
        else 
        {
          int color = m_is_new ? 1 : MPI_UNDEFINED;
          MPI_Comm_split(MPI_COMM_WORLD, color, 0, &m_comm_new);
        }
      } else
      {
        int comm_size_new;
        MPI_Comm_size(comm_new, &comm_size_new);
        assertAlways(comm_size_new == to_procs, "provided communicator has size different than number of procs to split mesh into");
        m_comm_new = comm_new;
      }
    }

    apf::Mesh2* m_mesh;
    gmi_model* m_model;
    MPI_Comm m_comm_new;      // the MPI_Comm over the given number of procs
    int m_myrank_world;
    int m_comm_size_world;
    bool m_is_new;

};

}

#endif