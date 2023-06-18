#ifndef HEAT_CALC_MESH_DIRICHLET_UPDATE_MAP
#define HEAT_CALC_MESH_DIRICHLET_UPDATE_MAP

#include <apfShape.h>
#include <vector>
#include "ProjectDefs.h"
#include "apf_data.h"
#include "mesh/apfMDSField.h"
#include "mesh/mesh_input.h"
#include "utils/mpi_utils.h"

namespace Mesh {

class MeshCG;

class DirichletUpdateMap
{
  public:
    DirichletUpdateMap(MeshCG* mesh, MPI_Comm comm);

    struct ArrayNode
    {
      int block;
      int el;
      int localnode;
    };

    int getNumLocalSrcNodes() const { return m_local_src_nodes.size(); }

    // returns the src nodes and populates dest_nodes with the destinations
    ArrayNode getLocalNodes(int idx, std::vector<ArrayNode>& dest_nodes)
    {
      dest_nodes.clear();
      for (int i=m_local_dest_node_ptrs[idx]; i < m_local_dest_node_ptrs[idx+1]; ++i)
        dest_nodes.push_back(m_local_dest_nodes[i]);

      return m_local_src_nodes[idx];
    }

    int getSendCount(int rank) const { return m_send_nodes[rank].size(); }

    int getRecvCount(int rank) const { return m_recv_node_ptrs[rank].size() - 1; }

    const std::vector<ArrayNode>& getSendDofs(int rank) const { return m_send_nodes[rank]; }

    // gives all the nodes that a value from the given rank and receive buffer
    // idx should be written to
    void getRecvDofs(int rank, int idx, std::vector<ArrayNode>& dofs) const
    {
      dofs.clear();
      for (int i=m_recv_node_ptrs[rank][idx]; i < m_recv_node_ptrs[rank][idx + 1]; ++i)
        dofs.push_back(m_recv_nodes[rank][i]);
    }

  private:
    static constexpr int NO_ENTITY = std::numeric_limits<int>::max();
    struct SurfaceModelEntity
    {
      apf::MeshEntity* vert = nullptr;
      int rank = 0;
      int model_entity_dim = NO_ENTITY;
      int model_entity_tag = NO_ENTITY;
    };

    struct MeshEntityStatus
    {
      apf::MeshEntity* e;
      bool status;
    };

    using ModelEntityField = fast_field::ApfMDSField<SurfaceModelEntity>;

    std::shared_ptr<ModelEntityField> collectMinDirichletSurfaces();

    SurfaceModelEntity getMinDirichletSurface(apf::MeshEntity* e);

    bool hasLocalDirichletSurface(apf::MeshEntity* e);

    bool isSurfaceDirichlet(const ModelEntitySpec& me_spec);

    void createSendAndRecvLists(std::shared_ptr<ModelEntityField> min_dirichlet_surface);

    //void unpackRecvEntities(ParallelExchange<apf::MeshEntity*>& exchanger);

    void getLocalDirichletUpdateMap();

    ArrayNode getArrayNode(apf::MeshEntity* vert);

    void getArrayNodes(apf::MeshEntity* vert, std::vector<ArrayNode>& nodes);

    ArrayNode getArrayNode(apf::MeshEntity* el, apf::MeshEntity* vert);


    MeshCG* m_mesh;
    const ApfData& m_apf_data;
    std::vector<MeshEntityGroupSpec> m_bc_specs;
    MPI_Comm m_comm;

    std::vector<ArrayNode> m_local_src_nodes;
    std::vector<ArrayNode> m_local_dest_nodes;
    std::vector<int> m_local_dest_node_ptrs;

    std::vector<std::vector<ArrayNode>> m_send_nodes;
    std::vector<int> m_recv_counts;
    std::vector<std::vector<ArrayNode>> m_recv_nodes;
    std::vector<std::vector<int>> m_recv_node_ptrs;



  friend bool operator<(const DirichletUpdateMap::SurfaceModelEntity& lhs, const DirichletUpdateMap::SurfaceModelEntity& rhs);


};

inline bool operator<(const DirichletUpdateMap::SurfaceModelEntity& lhs, const DirichletUpdateMap::SurfaceModelEntity& rhs)
{
  if (lhs.model_entity_dim != rhs.model_entity_dim)
  {
    return lhs.model_entity_dim < rhs.model_entity_dim;
  } else if (lhs.model_entity_tag != rhs.model_entity_tag)
  {
    return lhs.model_entity_tag < rhs.model_entity_tag;
  } else
  {
    return lhs.rank < rhs.rank;
  }
}

}

#endif