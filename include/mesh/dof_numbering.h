#ifndef A2_H
#define A2_H

#include <iostream>
#include <deque>
#include <queue>

#include "ProjectDefs.h"
#include "apf.h"
#include "apfNumbering.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "pumi.h"

#include "mesh/apfMDSField.h"
#include "mpi.h"

namespace Mesh{

// number dofs using adjacency-based algorithm
// in global mode, dof numbers are globally unique, and dirichlet nodes are all given
// the same value, the max value for int
class AdjacencyNumberer
{
#ifdef MESH_USE_MDS_NUMBERING
  using NumberingType = apf::ApfMDSNumbering;
#else
  using NumberingType = apf::Numbering;
#endif
  public:
    AdjacencyNumberer(apf::Mesh2* m, NumberingType* dof_nums,
                      NumberingType* el_nums, NumberingType* is_dirichlet,
                      bool global_numbering = false
                     ) :
      m_local(m),
      m_dof_nums(dof_nums),
      m_el_nums(el_nums),
      m_is_dirichlet(is_dirichlet),
      m_number_owned_nodes_only(global_numbering),
      m_number_dirichlet_constant(global_numbering)
    {
      m_ncomp = apf::countComponents(dof_nums);
      numberElements(m_local);
      countNodes(m_local, m_is_dirichlet, m_num_nodes, m_num_dirichlet);
      m_num_local_nodes = m_num_nodes;

      if (global_numbering)
      {
        m_local_owned_to_global_offset = 0;
        MPI_Exscan(&m_num_nodes, &m_local_owned_to_global_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        m_num_nodes += m_local_owned_to_global_offset;
      }
#ifdef MESH_USE_MDS_NUMBERING
      m_node_status = apf::createNumberingMDS(m_local, "dof_status",
                                              apf::getShape(dof_nums), 1);
#else
      m_node_status = apf::createNumbering(m_local, "dof_status",
                                           apf::getShape(dof_nums), 1);
#endif
      numberdofs(m_ncomp*(m_num_nodes + m_num_dirichlet), m_ncomp);

      if (global_numbering)
        apf::synchronize(dof_nums);
    }

    ~AdjacencyNumberer()
    {
      apf::destroyNumbering(m_node_status);
    }

    void reorder();

    int getNumDofs() const {return m_num_local_nodes * apf::countComponents(m_dof_nums);}

    int getNumTotalDofs() const {return (m_num_local_nodes + m_num_dirichlet) * apf::countComponents(m_dof_nums);}

    int getLocalOwnedToGlobalOffset() const { return m_local_owned_to_global_offset; }

    // returns vector of mesh elements in numbered order
    std::vector<apf::MeshEntity*> getElements();

  private:
    bool hasNode(apf::Mesh2* m, apf::MeshEntity* e);
    int nodeCount(apf::Mesh2* m_local, apf::MeshEntity* e);
    bool shouldNumber(NumberingType* m_is_dirichlet, apf::MeshEntity* e,
                      int node, int component);
    void addQueues(std::queue<apf::MeshEntity*> & q1, std::queue<apf::MeshEntity*> & q2);
    apf::MeshEntity* getStartEntity(apf::Mesh2* m);
    void printType(apf::Mesh* m, apf::MeshEntity* e);
    void countNodes(apf::Mesh2* m, NumberingType* is_dirichlet, int& num_nodes, int& num_dirichlet);
    void numberdofs(int ndof, int comp);
    bool isQueued(apf::MeshEntity* e, int node);
    bool isLabeled(apf::MeshEntity* e, int node);

    void numberElements(apf::Mesh2* m);
    void printElNumbers(apf::Mesh2* m);
    void labelNode(apf::MeshEntity* e, int node);
    void setQueued(apf::MeshEntity* e);

    enum
    {
      LABELED=0,
      QUEUED,
      UNSEEN,
    };


    apf::Mesh2* m_local;
    NumberingType* m_dof_nums;
    NumberingType* m_el_nums;
    NumberingType* m_is_dirichlet;
    NumberingType* m_node_status;
    bool m_number_owned_nodes_only = false;
    bool m_number_dirichlet_constant = false;
    int m_num_local_nodes = 0;
    int m_num_nodes = 0;
    int m_local_owned_to_global_offset = 0;
    int m_num_dirichlet = 0;
    int m_node_label = 0;
    int m_dirichlet_label = 0;
    int m_ncomp = 0;
};

} // namespace
#endif
