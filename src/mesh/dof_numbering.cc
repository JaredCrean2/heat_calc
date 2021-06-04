#include <climits>
#include <cassert>
#include <limits>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <deque>
#include <queue>

#include "mesh/dof_numbering.h"

namespace Mesh {

bool AdjacencyNumberer::hasNode(apf::Mesh2* m_local, apf::MeshEntity* e)
{
  return m_local->getShape()->countNodesOn(m_local->getType(e)) > 0;
}

int AdjacencyNumberer::nodeCount(apf::Mesh2* m_local, apf::MeshEntity* e)
{
  return m_local->getShape()->countNodesOn(m_local->getType(e));
}

// determine if a dof should be numbered or not
// Takes into account whether m_is_dirichlet is NULL or not
bool AdjacencyNumberer::shouldNumber(apf::Numbering* m_is_dirichlet, apf::MeshEntity* e,
                  int node, int component)
{
  bool ret = true;
  // number if dof is free or loaded
  if (m_is_dirichlet)
    ret = (!apf::getNumber(m_is_dirichlet, e, node, component));

  return ret;

}

// adds p2 to the end of q1
void AdjacencyNumberer::addQueues(std::queue<apf::MeshEntity*> & q1,
                                  std::queue<apf::MeshEntity*> & q2)
{
  int sizeq2 = q2.size();
  for ( int i = 0; i < sizeq2; ++i)
  {
    q1.push(q2.front());
    q2.pop();
  }
}

// get starting entity for node reordering
// search for node vertex classified on geometric vertex that is closest to 
// the point (x,y) and has minimum connectivity
//
apf::MeshEntity* AdjacencyNumberer::getStartEntity(apf::Mesh2* m_local)
{
  apf::MeshEntity* e_min; // minimum degree meshentity
  int degree_min = std::numeric_limits<int>::max();
  apf::MeshEntity* e_i;   // current meshentity

  apf::MeshIterator* it = m_local->begin(0); // iterator over verticies
  e_min = NULL;

  while ( (e_i = m_local->iterate(it)) )
  {
    int degree = m_local->countUpward(e_i);
    if (degree < degree_min)
    {
      e_min = e_i;
      degree_min = degree;
    }

  }

  assert(e_min);  // it shouldn't be possible to not find a an entity
  return e_min;
}


void AdjacencyNumberer::printType(apf::Mesh* m_local, apf::MeshEntity* e)
{
  int type_enum = m_local->getType(e);
  // check type, print appropriate message
    
  switch(type_enum)
    {
        case apf::Mesh::TET :
          std::cout << " tetrahedron" << std::endl;
          break;
        case apf::Mesh::HEX :
          std::cout << " Hexahedron" << std::endl;
           break;
        case apf::Mesh::PRISM :
           std::cout << " prism" << std::endl;
           break;
        case apf::Mesh::PYRAMID :
          std::cout << " pyramid" << std::endl;
          break;
        case apf::Mesh::QUAD :
          std::cout << " quad" << std::endl;
          break;
        case apf::Mesh::TRIANGLE :
          std::cout << " triangle" << std::endl;
          break;
        case apf::Mesh::VERTEX :
          std::cout << "vertex" << std::endl;
          break;
        case apf::Mesh::EDGE :
          std::cout << "edge" << std::endl;
          break;
        default:
             std::cout << " has no matching type" << std::endl;
    }
}

void AdjacencyNumberer::countNodes(apf::Mesh2* m_local, apf::Numbering* is_dirichlet, int& num_nodes, int& num_dirichlet)
{
  apf::FieldShape* fshape = apf::getShape(is_dirichlet);

  num_nodes = 0;
  num_dirichlet = 0;
  for (int dim=0; dim <= m_local->getDimension(); ++dim)
  {
    apf::MeshIterator* it = m_local->begin(dim);
    apf::MeshEntity* e;
    while ( (e = m_local->iterate(it)) )
    {
      int type = m_local->getType(e);
      for (int node=0; node < fshape->countNodesOn(type); ++node)
        if (!apf::getNumber(is_dirichlet, e, node, 0))
            num_nodes += 1;
        else
            num_dirichlet += 1;
    }
    m_local->end(it);
  }
}

// initially number all dofs with number greater than number of nodes, to show
// they have not received final number yet

void AdjacencyNumberer::numberdofs(int ndof, int comp)
{
//  apf::FieldShape* fieldshape = m_local->getShape();
  apf::MeshIterator* it;
  apf::MeshEntity* e;
  int numNodes_typei;

  for (int i = 0; i < 4; ++i) // loop over entity types
  {
    it = m_local->begin(i);
    e = m_local->deref(it);
    numNodes_typei = nodeCount(m_local,e);
    it = m_local->begin(i);

    if (numNodes_typei)  // if there are any dofs on this type of entity
      while ( (e = m_local->iterate(it)) )
        for ( int j = 0; j < numNodes_typei; ++j)
        {
          for ( int c = 0; c < comp; ++c)
          {
            apf::number(m_dof_nums, e, j, c, ndof); // label current node
          }
          apf::number(m_node_status, e, j, 0, UNSEEN);
        }
  }
}


bool AdjacencyNumberer::isQueued(apf::MeshEntity* e, int node)
{
  return apf::getNumber(m_node_status, e, node, 0) == QUEUED;
}

bool AdjacencyNumberer::isLabeled(apf::MeshEntity* e, int node)
{
  return apf::getNumber(m_node_status, e, node, 0) == LABELED;
}



// number elements with numbers greater than number of elements to show they 
// have not yet received their final number
void AdjacencyNumberer::numberElements(apf::Mesh2* m_local)
{

  if (m_el_nums)
  {

    const int numEl = m_local->count(m_local->getDimension());  // counts the number of elements
    apf::MeshIterator* it = m_local->begin(m_local->getDimension());
    apf::MeshEntity* e;
    int k = numEl + 1;

    while ( (e = m_local->iterate(it) ) )
    {
  //    std::cout << "labelling element " << k - numEl << " as " << k << std::endl;
      apf::number(m_el_nums, e, 0, 0, k);
      ++k;
    }
  }

}

void AdjacencyNumberer::printElNumbers(apf::Mesh2* m_local)
{
  apf::MeshIterator* it = m_local->begin(2);
  apf::MeshEntity* e;
  int i = 1;
  while ((e = m_local->iterate(it) ) )
  {
    int num = apf::getNumber(m_el_nums, e, 0 , 0);
    std::cout << "element " << i << "1 number = " << num << std::endl;
    ++i;
  }
}

void AdjacencyNumberer::labelNode(apf::MeshEntity* e, int node)
{
  for ( int c = 0; c < m_ncomp; ++c) // loop over dof of the node
  {
    if (shouldNumber(m_is_dirichlet, e, node, c))
      apf::number(m_dof_nums, e, node, c, --m_node_label);
    else
      apf::number(m_dof_nums, e, node, c, --m_dirichlet_label);
  }

  apf::number(m_node_status, e, node, 0, LABELED);
}

void AdjacencyNumberer::setQueued(apf::MeshEntity* e)
{
  for (int i=0; i < nodeCount(m_local, e); ++i)
    apf::number(m_node_status, e, i, 0, QUEUED);
}

// reorder mesh nodes and elements
// multiple dof on the same node are labelled sequentially
void AdjacencyNumberer::reorder()
{
// TODO: move m_is_dirichlet checks out one loop level because 
//       it is node status now, not dof status

  const int dim = m_local->getDimension();
  const int numEl = m_local->count(m_local->getDimension());  // counts the number of elements

  // initially number in range numnodes + 1 to 2*numnodes
  //numberElements(m_local, m_el_nums, numEl);
  //numberdofs(m_local, m_dof_nums, ndof, comp);


  // create queues
  std::queue <apf::MeshEntity*> que1;  // main queue
  std::queue <apf::MeshEntity*> tmpQue;  // temporary queue for vertices

  apf::MeshEntity* e;
  e = getStartEntity(m_local); // get starting node

  //int nodeLabel_i = ndof;
  m_node_label = m_num_nodes * m_ncomp;
  m_dirichlet_label = (m_num_nodes + m_num_dirichlet) * m_ncomp;
  int elementLabel_i = numEl;
  int numNodes_i;

  // queue initial entity
  que1.push(e);
 
  while (que1.size() > 0 )
  {
    e = que1.front(); // get next element in queue
    que1.pop();  // remove that element from the que

    // label all nodes on entity e
    numNodes_i = nodeCount(m_local,e);
    for ( int i = 0; i < numNodes_i; ++i)
      labelNode(e, i);

    // if e is a vertex, find adjacencies, look for unlabeled nodes
    if (m_local->getType(e) == apf::Mesh::VERTEX)
    {
      int numEdges_w = m_local->countUpward(e);
      for (int i = 0; i < numEdges_w; ++i)  // loop over edges
      {
        apf::MeshEntity* edge_i = m_local->getUpward(e, i);
        int numFaces_i = m_local->countUpward(edge_i);

        // queue faces adjacent to the edge if needed
        for (int j = 0; j < numFaces_i; ++j)  // loop over faces
        {
          apf::MeshEntity* face_j = m_local->getUpward(edge_i, j);

          if (m_el_nums && dim == 2)
          {
            // label face (element) if it is not yet labelled
            int faceNum_j = apf::getNumber(m_el_nums, face_j, 0, 0);
            if (faceNum_j > numEl) // if face not labelled
            {
              elementLabel_i -= 1;
              apf::number(m_el_nums, face_j, 0, 0, elementLabel_i);
            }
          }

          // add face to queue if it hasn't been labelled yet
          //std::cout << "queing faces" << std::endl;
          if ( hasNode(m_local, face_j) )
          {
            for (int node=0; node < nodeCount(m_local, face_j); ++node)
              if (!isLabeled(face_j, node) && !isQueued(face_j, node))
              {
                que1.push(face_j);
                setQueued(face_j);
              }
          }

          // get regions adjacent to the face
          if (dim == 3)
          {
            int numRegions_i = m_local->countUpward(face_j);
            for (int k = 0; k < numRegions_i; ++k)
            {
              apf::MeshEntity* region_k = m_local->getUpward(face_j, k);

              // label element
              if (m_el_nums)
              {
                int regionNum_k = apf::getNumber(m_el_nums, region_k, 0, 0);
                if (regionNum_k >= numEl)  // face not labelled
                {
                  elementLabel_i -= 1;
                  apf::number(m_el_nums, region_k, 0, 0, elementLabel_i);
                }
              }

              // add region to queue if it hasn't been labelled yet
              //std::cout << "queing regions" << std::endl;
              if ( hasNode(m_local, region_k) )
              {
                for (int node=0; node < nodeCount(m_local, region_k); ++node)
                  if (!isLabeled(region_k, node) && !isQueued(region_k, node))
                  {
                    que1.push(region_k);
                    setQueued(region_k);
                  }
              }
            }  // end region loop
          }
        }  // end face loop

        // look at other vertex on edge
        apf::MeshEntity* otherVertex = apf::getEdgeVertOppositeVert(m_local, edge_i, e);
        bool labeled = isLabeled(otherVertex, 0);
        bool queued  = isQueued(otherVertex, 0);

        if (hasNode(m_local, edge_i))
        {
          bool edgeNotLabeled = (!isLabeled(edge_i, 0) && !isQueued(edge_i, 0)); // edge not labeled nor in que


          if ((labeled || queued) && edgeNotLabeled)
          {
            for (int j = 0; j < nodeCount(m_local, edge_i); ++j)
              labelNode(edge_i, j);
          } else // add edge to queue first, othervertex second (via list)
          {
          
            if (edgeNotLabeled) 
            {
              que1.push(edge_i);
              setQueued(edge_i);
            }

            if ( (!labeled) && (!queued) )
            {
              tmpQue.push(otherVertex);
              setQueued(otherVertex);
            }
          }          

        } else  // if edge does not have node, deal with otherVertex
        {
          if ( (!labeled) && (!queued))  // if otherVertex not labelled
          {
            tmpQue.push(otherVertex);
            setQueued(otherVertex);
          }
        }

      }  // end edge loop

      // copy tmpQue into que1, empty tmpQue
      addQueues(que1, tmpQue);


      }  // end if (vertex)

      /*
      // either the algorithm is broken or there are disjoint sets of mesh
      // entities.  Assume the second and find a new starting entity
      if (nodeLabel_i != 0 && que1.size() == 0)
      {
        std::cout << "getting new starting entity" << std::endl;
        que1.push(getStartEntity(m_local));
      }
      */
    } // end while loop over queue


  if (m_el_nums)
  {
    if (elementLabel_i != 0)
    {
      std::cerr << "Warning: element numbering not sane" << std::endl;
      std::cerr << "final elementLabel_i = " << elementLabel_i << std::endl;
    } else
    {
      //std::cout << "element reordering is sane" << std::endl;
    }
  }

  if (m_node_label != 0)
  {
    std::cerr << "Warning: node numbering not sane" << std::endl;
    std::cerr << "final node label = " << m_node_label << std::endl;
  } else
  {
//    std::cout << "node reordering is sane" << std::endl;
  }

  if (m_dirichlet_label != m_ncomp * m_num_nodes)
  {
    std::cerr << "Warning: dirichlet node numbering not sane" << std::endl;
    std::cerr << "final dirichlet label = " << m_dirichlet_label << std::endl;

  }

  assert(m_node_label == 0);
  assert(m_dirichlet_label == m_ncomp * m_num_nodes);
}

std::vector<apf::MeshEntity*> AdjacencyNumberer::getElements()
{
  std::vector<apf::MeshEntity*> elements(m_local->count(3));
  apf::MeshIterator* it = m_local->begin(3);
  apf::MeshEntity* e;
  while ( (e = m_local->iterate(it)) )
    elements[apf::getNumber(m_el_nums, e, 0, 0)] = e;

  return elements;
}

} // namespace
