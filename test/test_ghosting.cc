#include "gtest/gtest.h"
#include <apfVector.h>
#include "mesh/ghosting.h"
#include "mesh/mesh_generator.h"
#include "mesh/mesh_input.h"
#include "mpi_utils.h"

namespace {

apf::Vector3 computeCentroid(apf::Mesh* mesh, apf::MeshEntity* e)
{
  if (mesh->getType(e) == apf::Mesh::Type::VERTEX)
  {
    apf::Vector3 pt;
    mesh->getPoint(e, 0, pt);
    return pt;
  } else
  {
    apf::Vector3 centroid;
    apf::Downward down;
    int ndown = mesh->getDownward(e, 0, down);
    for (int i=0; i < ndown; ++i)
    {
      centroid += computeCentroid(mesh, down[i]);
    }

    return centroid / ndown;
  }
}

apf::MeshEntity* get_closest_entity(apf::Mesh2* mesh, int dim, const apf::Vector3& pt, double tol)
{
  apf::MeshIterator* it = mesh->begin(dim);
  apf::MeshEntity* e;

  while ( (e = mesh->iterate(it)) )
  {
    apf::Vector3 centroid = computeCentroid(mesh, e);
    if (std::abs(pt.x() - centroid.x()) < tol &&
        std::abs(pt.y() - centroid.y()) < tol &&
        std::abs(pt.z() - centroid.z()) < tol)
    {
      mesh->end(it);
      return e;
    }
  }

  mesh->end(it);
  return nullptr;
}

void getGhostsAndSelf(apf::Mesh* mesh, apf::MeshEntity* e, apf::Copies& ghosts)
{
  int myrank = commRank(MPI_COMM_WORLD);
  ghosts.clear();
  mesh->getGhosts(e, ghosts);
  EXPECT_EQ(ghosts.count(myrank), 0);
  ghosts[myrank] = e;
}

void check_ghosting_symmetric(apf::Mesh* mesh)
{
  int myrank = commRank(MPI_COMM_WORLD);
  for (int dim=0; dim <= 3; ++dim)
  {
    PCU_Comm_Begin();
    apf::MeshIterator* it = mesh->begin(dim);
    apf::MeshEntity* e;
    apf::Copies ghosts;
    while ( (e = mesh->iterate(it)) )
    {
      if (mesh->isGhost(e) || mesh->isGhosted(e))
      {
        getGhostsAndSelf(mesh, e, ghosts);

        for (auto& p : ghosts)
        {
          int dest = p.first;
          if (dest != myrank)
          {
            PCU_Comm_Pack(dest, int(ghosts.size()));
            for (auto& p2 : ghosts)
            {
              PCU_Comm_Pack(dest, p2.first);
              PCU_Comm_Pack(dest, p2.second);
            }
          }
        }
      }
    }
    mesh->end(it);

    PCU_Comm_Send();

    apf::Copies ghosts_received;
    while (PCU_Comm_Listen())
    {
      while (!PCU_Comm_Unpacked())
      {
        int nvals = PCU_Comm_Unpack<int>();
        ghosts_received.clear();
        for (int i=0; i < nvals; ++i)
        {
          int rank = PCU_Comm_Unpack<int>();
          auto entity_recv = PCU_Comm_Unpack<apf::MeshEntity*>();
          ghosts_received.insert(std::make_pair(rank, entity_recv));
        }

        e = ghosts_received.at(myrank);
        getGhostsAndSelf(mesh, e, ghosts);

        EXPECT_EQ(ghosts.size(), ghosts_received.size());
        EXPECT_EQ(ghosts, ghosts_received);
      }
    }
  }

}
}

// Tests:
//   - single proc: verify it doesn't change anything
//   - 2 procs, 2x2 (need to know how the mesh is partitioned)
//     - check number of entities
//     - check coordinates of verts
//     - check ghosts are symmetric
//     - check number of ghosts for each vert is correct
//   - 4 procs, 4x4
//     - same as above

TEST(Ghosting, SingleProc)
{
  if (commSize(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 2, 0, 2, 0, 1, 2, 2, 1);
  apf::Mesh2* mesh = Mesh::make_parallel_mesh(spec, commSize(MPI_COMM_WORLD), Mesh::identity);

  Mesh::createGhosting(mesh, MPI_COMM_WORLD);

  EXPECT_EQ(mesh->count(0), 18);
  EXPECT_EQ(mesh->count(1), 2*12 + 9);
  EXPECT_EQ(mesh->count(2), 2*4 + 12);
  EXPECT_EQ(mesh->count(3), 4);
}

TEST(Ghosting, 2Procs)
{
  if (commSize(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 2, 0, 2, 0, 1, 2, 2, 1);
  apf::Mesh2* mesh = Mesh::make_parallel_mesh(spec, commSize(MPI_COMM_WORLD), Mesh::identity);

  Mesh::createGhosting(mesh, MPI_COMM_WORLD);

  EXPECT_EQ(mesh->count(0), 18);
  EXPECT_EQ(mesh->count(1), 2*12 + 9);
  EXPECT_EQ(mesh->count(2), 2*4 + 12);
  EXPECT_EQ(mesh->count(3), 4);

  std::vector<apf::Vector3> edge_centroids =  { // parallel to y axis
                                                {0, 0.5, 0}, {0, 1.5, 0},
                                                {1, 0.5, 0}, {1, 1.5, 0},
                                                {2, 0.5, 0}, {2, 1.5, 0},
                                                // parallel to x axis
                                                {0.5, 0, 0}, {1.5, 0, 0},
                                                {0.5, 1, 0}, {1.5, 1, 0},
                                                {0.5, 2, 0}, {1.5, 2, 0}
                                              };
                                              

  std::vector<apf::Vector3> face_centroids = { // xy plane, bottom
                                               {0.5, 0.5, 0}, {1.5, 0.5, 0},
                                               {0.5, 1.5, 0}, {1.5, 1.5, 0},
                                               // xy plane, top
                                               {0.5, 0.5, 1}, {1.5, 0.5, 1},
                                               {0.5, 1.5, 1}, {1.5, 1.5, 1},

                                               {0.5, 0, 0.5}, {1.5, 0, 0.5},
                                               {0.5, 1, 0.5}, {1.5, 1, 0.5},
                                               {0.5, 2, 0.5}, {1.5, 2, 0.5},

                                               {0, 0.5, 0.5}, {0, 1.5, 0.5},
                                               {1, 0.5, 0.5}, {1, 1.5, 0.5},
                                               {2, 0.5, 0.5}, {2, 1.5, .5}
                                             };           


  std::vector<apf::Vector3> el_centroids = { {0.5, 0.5, 0.5}, {1.5, 0.5, 0.5},
                                             {0.5, 1.5, 0.5}, {1.5, 1.5, 0.5}
                                           };                                        


  // each triplet is (isShared, isGhost, isGhosted)
  std::vector<std::array<bool, 3>> vert_desc, edge_desc, face_desc, el_desc;
  if (commRank(MPI_COMM_WORLD) == 0)
  {
    vert_desc = { {false, false, true}, {false, false, true}, {false, false, true},
                  {true, false, false}, {true, false, false}, {true, false, false},
                  {false, true, false}, {false, true, false}, {false, true, false},
                };

    edge_desc = { {false, false, true}, {false, true, false},
                  {false, false, true}, {false, true, false},
                  {false, false, true}, {false, true, false},
                  
                  {false, false, true}, {false, false, true},
                  {true, false, false}, {true, false, false},
                  {false, true, false}, {false, true, false}
                };

    face_desc = { {false, false, true}, {false, false, true},
                  {false, true, false}, {false, true, false},

                  {false, false, true}, {false, false, true},
                  {false, true, false}, {false, true, false},

                  {false, false, true}, {false, false, true},
                  {true, false, false}, {true, false, false},
                  {false, true, false}, {false, true, false},

                  {false, false, true}, {false, true, false},
                  {false, false, true}, {false, true, false},
                  {false, false, true}, {false, true, false}
                };

    el_desc = { {false, false, true}, {false, false, true},
                {false, true, false}, {false, true, false}
              };                

  } else
  {
    vert_desc = { {false, true, false}, {false, true, false}, {false, true, false},
                  {true, false, false}, {true, false, false}, {true, false, false},
                  {false, false, true}, {false, false, true}, {false, false, true},
                }; 

    edge_desc = { {false, true, false}, {false, false, true},
                  {false, true, false}, {false, false, true},
                  {false, true, false}, {false, false, true},
                  
                  {false, true, false}, {false, true, false},
                  {true, false, false}, {true, false, false},
                  {false, false, true}, {false, false, true}
                };

    face_desc = { {false, true, false}, {false, true, false},
                  {false, false, true}, {false, false, true},

                  {false, true, false}, {false, true, false},
                  {false, false, true}, {false, false, true},

                  {false, true, false}, {false, true, false},
                  {true, false, false}, {true, false, false},
                  {false, false, true}, {false, false, true},


                  {false, true, false}, {false, false, true},
                  {false, true, false}, {false, false, true},
                  {false, true, false}, {false, false, true}
                };

    el_desc = { {false, true, false}, {false, true, false},
                {false, false, true}, {false, false, true},
              };                 
  }

  for (int k=0; k < spec.nz+1; ++k)
  {
    int idx = 0; // the z=0 verts follow the same pattern as the z=1 verts
    for (int j=0; j < spec.ny+1; ++j)
        for (int i=0; i < spec.nx+1; ++i)
        {
          double x = i, y = j, z = k;
          apf::MeshEntity* vert = get_closest_entity(mesh, 0, {x, y, z}, 1e-12);
          EXPECT_TRUE(vert);
          EXPECT_EQ(mesh->isShared(vert),  vert_desc[idx][0]);
          EXPECT_EQ(mesh->isGhost(vert),   vert_desc[idx][1]);
          EXPECT_EQ(mesh->isGhosted(vert), vert_desc[idx][2]);
          idx++;
        }
  } 

  for (size_t i=0; i < edge_centroids.size(); ++i)
  {
    apf::MeshEntity* edge = get_closest_entity(mesh, 1, edge_centroids[i], 1e-12);

    EXPECT_TRUE(edge); 
    EXPECT_EQ(mesh->isShared(edge),  edge_desc[i][0]);
    EXPECT_EQ(mesh->isGhost(edge),   edge_desc[i][1]);
    EXPECT_EQ(mesh->isGhosted(edge), edge_desc[i][2]);      
  }

  for (size_t i=0; i < face_centroids.size(); ++i)
  {
    apf::MeshEntity* face = get_closest_entity(mesh, 2, face_centroids[i], 1e-12);
    
    EXPECT_TRUE(face); 
    EXPECT_EQ(mesh->isShared(face),  face_desc[i][0]);
    EXPECT_EQ(mesh->isGhost(face),   face_desc[i][1]);
    EXPECT_EQ(mesh->isGhosted(face), face_desc[i][2]);      
  } 

  for (size_t i=0; i < el_centroids.size(); ++i)
  {
    apf::MeshEntity* el = get_closest_entity(mesh, 3, el_centroids[i], 1e-12);
    
    EXPECT_TRUE(el); 
    EXPECT_EQ(mesh->isShared(el),  el_desc[i][0]);
    EXPECT_EQ(mesh->isGhost(el),   el_desc[i][1]);
    EXPECT_EQ(mesh->isGhosted(el), el_desc[i][2]);      
  }  

  check_ghosting_symmetric(mesh); 
}


TEST(Ghosting, 4Procs)
{
  if (commSize(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 4, 0, 4, 0, 1, 4, 4, 1);
  apf::Mesh2* mesh = Mesh::make_parallel_mesh(spec, commSize(MPI_COMM_WORLD), Mesh::identity);

  Mesh::createGhosting(mesh, MPI_COMM_WORLD);

  EXPECT_EQ(mesh->count(0), 32);
  EXPECT_EQ(mesh->count(1), 2*2*3*4 + 4*4);
  EXPECT_EQ(mesh->count(2), 3*3*2 + 2*3*4);
  EXPECT_EQ(mesh->count(3), 9);

  double x0, y0, z0;
  int myrank = commRank(MPI_COMM_WORLD);
  // each triplet is (isShared, isGhost, isGhosted)
  std::vector<std::array<bool, 3>> vert_desc, edge_desc, face_desc, el_desc;
  if (myrank == 0)
  {
    x0 = 0;
    y0 = 0;
    z0 = 0;
    vert_desc = { {false, false, false}, {false, false, true}, {true, false, false}, {false, true, false},
                  {false, false, true},  {false, false, true}, {true, false, true},  {false, true, false},
                  {true, false, false},  {true, false, true},  {true, false, false}, {false, true, false},
                  {false, true, false},  {false, true, false}, {false, true, false}, {false, true, false}
                };

    edge_desc = { {false, false, false}, {false, false, true}, {false, true, false},
                  {false, false, true},  {false, false, true}, {false, true, false},
                  {true, false, false},  {true, false, true},  {false, true, false},
                  {false, true, false},  {false, true, false}, {false, true, false},
                  
                  {false, false, false}, {false, false, true}, {false, true, false},
                  {false, false, true},  {false, false, true}, {false, true, false},
                  {true, false, false},  {true, false, true},  {false, true, false},
                  {false, true, false},  {false, true, false}, {false, true, false},

                  {false, false, false}, {false, false, true}, {true, false, false}, {false, true, false},
                  {false, false, true}, {false, false, true}, {true, false, true}, {false, true, false},
                  {true, false, false},  {true, false, true},  {true, false, false}, {false, true, false},
                  {false, true, false},  {false, true, false}, {false, true, false}, {false, true, false}
                };
                
  } else if (myrank == 1)
  {
    x0 = 0;
    y0 = 1;
    z0 = 0;

    vert_desc = { {false, true, false}, {false, true, false}, {false, true, false}, {false, true, false},
                  {true, false, false},  {true, false, true}, {true, false, false},  {false, true, false},
                  {false, false, true},  {false, false, true},  {true, false, true}, {false, true, false},
                  {false, false, false},  {false, false, true}, {true, false, false}, {false, true, false}
                };

    edge_desc = { {false, true, false},  {false, true, false}, {false, true, false},
                  {true, false, false},  {true, false, true},  {false, true, false},
                  {false, false, true},  {false, false, true}, {false, true, false},
                  {false, false, false}, {false, false, true}, {false, true, false},
                  
                  {false, true, false},  {false, false, true}, {false, false, false},
                  {false, true, false},  {false, false, true}, {false, false, true},
                  {false, true, false},  {true, false, true},  {true, false, false},
                  {false, true, false},  {false, true, false}, {false, true, false},

                  {false, true, false}, {false, true, false}, {false, true, false}, {false, true, false},
                  {true, false, false}, {true, false, true}, {true, false, false}, {false, true, false},
                  {false, false, true},  {false, false, true},  {true, false, true}, {false, true, false},
                  {false, false, false},  {false, false, true}, {true, false, false}, {false, true, false}
                };    
  } else if (myrank == 2)
  {
    x0 = 1;
    y0 = 0;
    z0 = 0;

    vert_desc = { {false, true, false}, {true, false, false}, {false, false, true}, {false, false, false},
                  {false, true, false},  {true, false, true}, {false, false, true},  {false, false, true},
                  {false, true, false},  {true, false, false},  {true, false, true}, {true, false, false},
                  {false, true, false},  {false, true, false}, {false, true, false}, {false, true, false}
                };

    edge_desc = { {false, true, false}, {false, false, true}, {false, false, false},
                  {false, true, false}, {false, false, true}, {false, false, true},
                  {false, true, false}, {true, false, true},  {true, false, false},
                  {false, true, false}, {false, true, false}, {false, true, false},
                  
                  {false, true, false}, {false, true, false}, {false, true, false},
                  {true, false, false},  {true, false, true}, {false, true, false},
                  {false, false, true},  {false, false, true}, {false, true, false},
                  {false, false, false},  {false, false, true}, {false, true, false},

                  {false, true, false}, {true, false, false}, {false, false, true}, {false, false, false},
                  {false, true, false}, {true, false, true}, {false, false, true}, {false, false, true},
                  {false, true, false},  {true, false, false},  {true, false, true}, {true, false, false},
                  {false, true, false},  {false, true, false}, {false, true, false}, {false, true, false}
                };    
  } else
  {
    x0 = 1;
    y0 = 1;
    z0 = 0;

    vert_desc = { {false, true, false}, {false, true, false}, {false, true, false}, {false, true, false},
                  {false, true, false},  {true, false, false}, {true, false, true},  {true, false, false},
                  {false, true, false},  {true, false, true},  {false, false, true}, {false, false, true},
                  {false, true, false},  {true, false, false}, {false, false, true}, {false, false, false}
                };

    edge_desc = { {false, true, false}, {false, true, false}, {false, true, false},
                  {false, true, false},  {true, false, true}, {true, false, false},
                  {false, true, false},  {false, false, true},  {false, false, true},
                  {false, true, false},  {false, false, true}, {false, false, false},
                  
                  {false, true, false}, {false, true, false}, {false, true, false},
                  {false, true, false},  {true, false, true}, {true, false, false},
                  {false, true, false},  {false, false, true}, {false, false, true},
                  {false, true, false},  {false, false, true}, {false, false, false},

                  {false, true, false}, {false, true, false}, {false, true, false}, {false, true, false},
                  {false, true, false}, {true, false, false}, {true, false, true},  {true, false, false},
                  {false, true, false},  {true, false, true},  {false, false, true}, {false, false, true},
                  {false, true, false},  {true, false, false}, {false, false, true}, {false, false, false}
                };
  }

  std::vector<apf::Vector3> edge_centroids = { // parallel to x axis 
                                               {0.5, 0, 0}, {1.5, 0, 0}, {2.5, 0, 0},
                                               {0.5, 1, 0}, {1.5, 1, 0}, {2.5, 1, 0},
                                               {0.5, 2, 0}, {1.5, 2, 0}, {2.5, 2, 0},
                                               {0.5, 3, 0}, {1.5, 3, 0}, {2.5, 3, 0},
                                               // parallel to y axis
                                               {0, 0.5, 0}, {0, 1.5, 0}, {0, 2.5, 0},
                                               {1, 0.5, 0}, {1, 1.5, 0}, {1, 2.5, 0},
                                               {2, 0.5, 0}, {2, 1.5, 0}, {2, 2.5, 0},
                                               {3, 0.5, 0}, {3, 1.5, 0}, {3, 2.5, 0},
                                               // parallel to z axis
                                               {0, 0, 0.5}, {1, 0, 0.5}, {2, 0, 0.5}, {3, 0, 0.5},
                                               {0, 1, 0.5}, {1, 1, 0.5}, {2, 1, 0.5}, {3, 1, 0.5},
                                               {0, 2, 0.5}, {1, 2, 0.5}, {2, 2, 0.5}, {3, 2, 0.5},
                                               {0, 3, 0.5}, {1, 3, 0.5}, {2, 3, 0.5}, {3, 3, 0.5},                                               
                                              };

  std::vector<apf::Vector3> face_centroids = {  // xy plane (bottom)
                                               {0.5, 0.5, 0}, {1.5, 0.5, 0}, {2.5, 0.5, 0},
                                               {0.5, 1.5, 0}, {1.5, 1.5, 0}, {2.5, 1.5, 0},
                                               {0.5, 2.5, 0}, {1.5, 2.5, 0}, {2.5, 2.5, 0},
                                               // xy plane (top)
                                               {0.5, 0.5, 1}, {1.5, 0.5, 1}, {2.5, 0.5, 1},
                                               {0.5, 1.5, 1}, {1.5, 1.5, 1}, {2.5, 1.5, 1},
                                               {0.5, 2.5, 1}, {1.5, 2.5, 1}, {2.5, 2.5, 1},
                                               // xz plane
                                               {0.5, 0, 0.5}, {1.5, 0, 0.5}, {2.5, 0, 0.5},
                                               {0.5, 1, 0.5}, {1.5, 1, 0.5}, {2.5, 1, 0.5},
                                               {0.5, 2, 0.5}, {1.5, 2, 0.5}, {2.5, 2, 0.5},
                                               {0.5, 3, 0.5}, {1.5, 3, 0.5}, {2.5, 3, 0.5},
                                               // yz plane
                                               {0, 0.5, 0.5}, {0, 1.5, 0.5}, {0, 2.5, 0.5},
                                               {1, 0.5, 0.5}, {1, 1.5, 0.5}, {1, 2.5, 0.5},
                                               {2, 0.5, 0.5}, {2, 1.5, 0.5}, {2, 2.5, 0.5},
                                               {3, 0.5, 0.5}, {3, 1.5, 0.5}, {3, 2.5, 0.5},                                              
                                             };

  std::vector<apf::Vector3> el_centroids = { {0.5, 0.5, 0.5}, {1.5, 0.5, 0.5}, {2.5, 0.5, 0.5},
                                             {0.5, 1.5, 0.5}, {1.5, 1.5, 0.5}, {2.5, 1.5, 0.5},
                                             {0.5, 2.5, 0.5}, {1.5, 2.5, 0.5}, {2.5, 2.5, 0.5},
                                           };

  for (auto& pt : edge_centroids)
    pt += apf::Vector3(x0, y0, z0);

  for (auto& pt : face_centroids)
    pt += apf::Vector3(x0, y0, z0);

  for (auto& pt : el_centroids)
    pt += apf::Vector3(x0, y0, z0);
   

  for (int k=0; k < 2; ++k)
  {
    int idx = 0; // the z=0 verts follow the same pattern as the z=1 verts
    for (int j=0; j < 4; ++j)
        for (int i=0; i < 4; ++i)
        {
          double x = i + x0, y = j + y0, z = k + z0;
          apf::MeshEntity* vert = get_closest_entity(mesh, 0, {x, y, z}, 1e-12);
          EXPECT_TRUE(vert);
          EXPECT_EQ(mesh->isShared(vert),  vert_desc[idx][0]);
          EXPECT_EQ(mesh->isGhost(vert),   vert_desc[idx][1]);
          EXPECT_EQ(mesh->isGhosted(vert), vert_desc[idx][2]);
          idx++;
        }
  } 

  for (size_t i=0; i < edge_centroids.size(); ++i)
  {
    apf::MeshEntity* edge = get_closest_entity(mesh, 1, edge_centroids[i], 1e-12);

    EXPECT_TRUE(edge); 
    EXPECT_EQ(mesh->isShared(edge),  edge_desc[i][0]);
    EXPECT_EQ(mesh->isGhost(edge),   edge_desc[i][1]);
    EXPECT_EQ(mesh->isGhosted(edge), edge_desc[i][2]);      
  }

  check_ghosting_symmetric(mesh);
}