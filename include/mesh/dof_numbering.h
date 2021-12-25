#ifndef A2_H
#define A2_H

#include <iostream>
#include <deque>
#include <queue>

#include "apf.h"
#include "apfNumbering.h"
#include "apfMesh2.h"
#include "apfShape.h"

namespace Mesh{

// number dofs using adjacency-based algorithm
class AdjacencyNumberer
{
  public:
    AdjacencyNumberer(apf::Mesh2* m, apf::Numbering* dof_nums,
                      apf::Numbering* el_nums, apf::Numbering* is_dirichlet
                     ) :
      m_local(m),
      m_dof_nums(dof_nums),
      m_el_nums(el_nums),
      m_is_dirichlet(is_dirichlet)
    {
      m_ncomp = apf::countComponents(dof_nums);
      numberElements(m_local);
      countNodes(m_local, m_is_dirichlet, m_num_nodes, m_num_dirichlet);
      m_node_status = apf::createNumbering(m_local, "dof_status",
                                           apf::getShape(dof_nums), 1);
      numberdofs(m_ncomp*(m_num_nodes + m_num_dirichlet), m_ncomp);
    }

    ~AdjacencyNumberer()
    {
      apf::destroyNumbering(m_node_status);
    }

    void reorder();

    int getNumDofs() const {return m_num_nodes * apf::countComponents(m_dof_nums);}

    int getNumTotalDofs() const {return (m_num_nodes + m_num_dirichlet) * apf::countComponents(m_dof_nums);}

    // returns vector of mesh elements in numbered order
    std::vector<apf::MeshEntity*> getElements();

  private:
    bool hasNode(apf::Mesh2* m, apf::MeshEntity* e);
    int nodeCount(apf::Mesh2* m_local, apf::MeshEntity* e);
    bool shouldNumber(apf::Numbering* m_is_dirichlet, apf::MeshEntity* e,
                      int node, int component);
    void addQueues(std::queue<apf::MeshEntity*> & q1, std::queue<apf::MeshEntity*> & q2);
    apf::MeshEntity* getStartEntity(apf::Mesh2* m);
    void printType(apf::Mesh* m, apf::MeshEntity* e);
    void countNodes(apf::Mesh2* m, apf::Numbering* is_dirichlet, int& num_nodes, int& num_dirichlet);
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
    apf::Numbering* m_dof_nums;
    apf::Numbering* m_el_nums;
    apf::Numbering* m_is_dirichlet;
    apf::Numbering* m_node_status;
    int m_num_nodes = 0;
    int m_num_dirichlet = 0;
    int m_node_label = 0;
    int m_dirichlet_label = 0;
    int m_ncomp = 0;
};

} // namespace
#endif
