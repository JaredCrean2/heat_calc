// generates structured meshes on a mapped cube

#ifndef MESH_GENERATOR_MULTI_H
#define MESH_GENERATOR_MULTI_H

#include <iostream>
#include "ProjectDefs.h"
#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "error_handling.h"
#include "gmi_null.h"
#include "mesh/mesh_input.h"
#include "mesh/mesh_generator.h"
#include "mesh/gmiDataStructure.h"

namespace Mesh
{


// creates a mesh that can be mapped to a series of cubes cube
// each block is defined by one MeshSpec, and the block share x-z planes
// face geometry classification:
//    0 = yminus
//    1 = xplus
//    2 = yplus
//    3 = xminus
//    4 = zminus
//    5 = zplus
template <typename T>
class MeshGeneratorMulti
{
  public:
    // func must be Point = func(const Point&)
    MeshGeneratorMulti(const std::vector<MeshSpec>& mesh_specs, T func) :
      m_func(func),
      m_meshspecs(mesh_specs),
      m_meshspec(nullptr),
      m_block_geometric_verts(8),
      m_block_geometric_edges(12),
      m_block_geometric_faces(6),
      m_block_geometric_entities{std::ref(m_block_geometric_verts), std::ref(m_block_geometric_edges), std::ref(m_block_geometric_faces)}
  {
    checkInput();
    int nx=0, ny=0, nz=0;
    for (auto& meshspec : m_meshspecs)
    {
      nx += meshspec.nx;
      ny += meshspec.ny;
      nz += meshspec.nz;
    }

    m_verts.resize(boost::extents[nx+1][ny+1][nz+1]);
  }

    //TODO: move function definitions to source file
    apf::Mesh2* generate()
    {
      gmi_register_null();
      m_gmi_topo = std::make_shared<mesh_gmi::GMITopo>();
      createGeometry();
      m_gmi_model = mesh_gmi::createGMITopo(m_gmi_topo);
      //m_gmi_model = gmi_load(".null");
      m_mesh = apf::makeEmptyMdsMesh(m_gmi_model, 3, false);

      createVertices();
      createEdges();
      createFaces();
      createElements();

      m_mesh->acceptChanges();
      m_mesh->verify();

      apf::writeASCIIVtkFiles("mesh_created", m_mesh);
     
      return m_mesh;
    }



  private:
    void createVertices()
    {
      for (unsigned int block=0; block < m_meshspecs.size(); ++block)
      {
        setBlock(block);
        int jstart = block == 0 ? 0 : 1;

        for (int i=0; i < m_meshspec->nx+1; ++i)
          for (int j=jstart; j < m_meshspec->ny+1; ++j)
            for (int k=0; k < m_meshspec->nz+1; ++k)
            {
              int j2 = j + m_yidx_offset;

              Point p = computeVertCoords(i, j, k);
              apf::Vector3 p2(p.x, p.y, p.z);
              auto me = getVertModelEntity(i, j, k);
              m_verts[i][j2][k] = m_mesh->createVert(me);
              m_mesh->setPoint(m_verts[i][j2][k], 0, p2);
            }
      }
    }

    apf::ModelEntity* getVertModelEntity(const int i, const int j, const int k)
    {
      int xmax = m_meshspec->nx;
      int ymax = m_meshspec->ny;
      int zmax = m_meshspec->nz;

      ModelEntitySpec me;
      // findModelEntity(dim, tag)
      // z = zmin vertices
      if (i == 0 && j == 0 && k == 0)
        me = ModelEntitySpec(0, m_block_geometric_verts[0]);
      else if (i == xmax && j == 0 && k == 0)
        me = ModelEntitySpec(0, m_block_geometric_verts[1]);
      else if (i == xmax && j == ymax && k == 0)
        me = ModelEntitySpec(0, m_block_geometric_verts[2]);
      else if (i == 0 && j == ymax && k == 0)
        me = ModelEntitySpec(0, m_block_geometric_verts[3]);
      // z = zmax vertices
      else if (i == 0 && j == 0 && k == zmax)
        me = ModelEntitySpec(0, m_block_geometric_verts[4]);
      else if (i == xmax && j == 0 && k == zmax)
        me = ModelEntitySpec(0, m_block_geometric_verts[5]);
      else if (i == xmax && j == ymax && k == zmax)
        me = ModelEntitySpec(0, m_block_geometric_verts[6]);
      else if (i == 0 && j == ymax && k == zmax)
        me = ModelEntitySpec(0, m_block_geometric_verts[7]);
      // edges on z = zmin quad
      else if (j == 0 && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[0]);
      else if (i == xmax && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[1]);
      else if (j == ymax && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[2]);
      else if (i == 0 && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[3]);
      // edges on z = zmax quad
      else if (j == 0 && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[4]);
      else if (i == xmax && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[5]);
      else if (j == ymax && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[6]);
      else if (i == 0 && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[7]);
      // edges parallel to z axis
      else if (i == 0 && j == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[8]);
      else if (i == xmax && j == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[9]);
      else if (i == xmax && j == ymax)
        me = ModelEntitySpec(1, m_block_geometric_edges[10]);
      else if (i == 0 && j == ymax)
        me = ModelEntitySpec(1, m_block_geometric_edges[11]);
      // faces
      else if (j == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[0]);
      else if (i == xmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[1]);
      else if (j == ymax)
        me = ModelEntitySpec(2, m_block_geometric_faces[2]);
      else if (i == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[3]);
      else if (k == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[4]);
      else if (k == zmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[5]);
      // interior
      else
        me = ModelEntitySpec(3, m_current_block);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    Point computeVertCoords(const int i, const int j, const int k)
    {
      double dx = (m_meshspec->xmax - m_meshspec->xmin)/m_meshspec->nx;
      double dy = (m_meshspec->ymax - m_meshspec->ymin)/m_meshspec->ny;
      double dz = (m_meshspec->zmax - m_meshspec->zmin)/m_meshspec->nz;

      Point p{i*dx + m_coords_min[0], j*dy + m_coords_min[1], k*dz + m_coords_min[2]};
      return m_func(p);
    }

    void createEdges()
    {
      for (unsigned int block=0; block < m_meshspecs.size(); ++block)
      {
        setBlock(block);
        int jstart = block == 0 ? 0 : 1;
        apf::MeshEntity* verts[2];
        apf::ModelEntity* me;
        for (int i=0; i < m_meshspec->nx; ++i)
          for (int j=jstart; j < m_meshspec->ny+1; ++j)
            for (int k=0; k < m_meshspec->nz+1; ++k)
            {
              // x direction
              int j2 = j + m_yidx_offset;
              verts[0] = m_verts[i][j2][k];
              verts[1] = m_verts[i+1][j2][k];
              me = getEdgeModelEntity_x(i, j, k);
              m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
            }

        for (int i=0; i < m_meshspec->nx+1; ++i)
          for (int j=0; j < m_meshspec->ny; ++j)
            for (int k=0; k < m_meshspec->nz+1; ++k)
            {
              // y direction
              int j2 = j + m_yidx_offset;
              verts[0] = m_verts[i][j2][k];
              verts[1] = m_verts[i][j2+1][k];
              me = getEdgeModelEntity_y(i, j, k);
              m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
            }

        for (int i=0; i < m_meshspec->nx+1; ++i)
          for (int j=jstart; j < m_meshspec->ny+1; ++j)
            for (int k=0; k < m_meshspec->nz; ++k)
            {
              // z direction
              int j2 = j + m_yidx_offset;
              verts[0] = m_verts[i][j2][k];
              verts[1] = m_verts[i][j2][k+1];
              me = getEdgeModelEntity_z(i, j, k);
              m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
            }
      }
    }

    // gets model entity for edges parallel to the x axis
    apf::ModelEntity* getEdgeModelEntity_x(const int i, const int j, const int k)
    {
      //int xmax = m_meshspec->nx + 1;
      int ymax = m_meshspec->ny;
      int zmax = m_meshspec->nz;

      ModelEntitySpec me;
      // geometric edges
      if (j == 0 && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[0]);
      else if (j == ymax && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[2]);
      else if (j == 0 && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[4]);
      else if (j == ymax && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[6]);
      // faces
      else if (j == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[0]);
      else if (k == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[4]);
      else if (j == ymax)
        me = ModelEntitySpec(2, m_block_geometric_faces[2]);
      else if (k == zmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[5]);
      // interior
      else
        me = ModelEntitySpec(3, m_current_block);


      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for edges parallel to the y axis
    apf::ModelEntity* getEdgeModelEntity_y(const int i, const int j, const int k)
    {
      int xmax = m_meshspec->nx;
      //int ymax = m_meshspec->ny;
      int zmax = m_meshspec->nz;

      ModelEntitySpec me;
      // geometric edges
      if (i == 0 && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[3]);
      else if (i == xmax && k == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[1]);
      else if (i == 0 && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[7]);
      else if (i == xmax && k == zmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[5]);
      // faces
      else if (i == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[3]);
      else if (k == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[4]);
      else if (i == xmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[1]);
      else if (k == zmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[5]);
      // interior
      else
        me = ModelEntitySpec(3, m_current_block);
 
      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for edges parallel to the z axis
    apf::ModelEntity* getEdgeModelEntity_z(const int i, const int j, const int k)
    {
      int xmax = m_meshspec->nx;
      int ymax = m_meshspec->ny;
      //int zmax = m_meshspec->nz;

      ModelEntitySpec me;
      // geometric edges
      if (j == 0 && i == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[8]);
      else if (j == ymax && i == 0)
        me = ModelEntitySpec(1, m_block_geometric_edges[11]);
      else if (j == 0 && i == xmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[9]);
      else if (j == ymax && i == xmax)
        me = ModelEntitySpec(1, m_block_geometric_edges[10]);
      // faces
      else if (j == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[0]);
      else if (i == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[3]);
      else if (j == ymax)
        me = ModelEntitySpec(2, m_block_geometric_faces[2]);
      else if (i == xmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[1]);
      // interior
      else
        me = ModelEntitySpec(3, m_current_block);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    void createFaces()
    {
      for (int block=0; block < m_meshspecs.size(); ++block)
      {
        setBlock(block);
        int jstart = block == 0 ? 0 : 1;
        apf::MeshEntity* verts[4];
        apf::ModelEntity* me;
        for (int i=0; i < m_meshspec->nx; ++i)
          for (int j=jstart; j < m_meshspec->ny+1; ++j)
            for (int k=0; k < m_meshspec->nz; ++k)
            {
              // xz plane
              int j2 = j + m_yidx_offset;
              verts[0] = m_verts[i][j2][k];
              verts[1] = m_verts[i+1][j2][k];
              verts[2] = m_verts[i+1][j2][k+1];
              verts[3] = m_verts[i][j2][k+1];
              me = getFaceModelEntity_xz(i, j, k);
              apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
            }

        for (int i=0; i < m_meshspec->nx+1; ++i)
          for (int j=0; j < m_meshspec->ny; ++j)
            for (int k=0; k < m_meshspec->nz; ++k)
            {
              // yz plane
              int j2 = j + m_yidx_offset;
              verts[0] = m_verts[i][j2][k];
              verts[1] = m_verts[i][j2][k+1];
              verts[2] = m_verts[i][j2+1][k+1];
              verts[3] = m_verts[i][j2+1][k];
              me = getFaceModelEntity_yz(i, j, k);
              apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
            }

        for (int i=0; i < m_meshspec->nx; ++i)
          for (int j=0; j < m_meshspec->ny; ++j)
            for (int k=0; k < m_meshspec->nz+1; ++k)
            {
              // xy plane
              int j2 = j + m_yidx_offset;
              verts[0] = m_verts[i][j2][k];
              verts[1] = m_verts[i+1][j2][k];
              verts[2] = m_verts[i+1][j2+1][k];
              verts[3] = m_verts[i][j2+1][k];
              me = getFaceModelEntity_xy(i, j, k);
              apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
            }
      }

    }

    // gets model entity for faces parallel to the xz plane
    apf::ModelEntity* getFaceModelEntity_xz(const int i, const int j, const int k)
    {
      int ymax = m_meshspec->ny;

      ModelEntitySpec me;
      if (j == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[0]);
      else if (j == ymax)
        me = ModelEntitySpec(2, m_block_geometric_faces[2]);
      else
        me = ModelEntitySpec(3, m_current_block);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for faces parallel to the yz plane
    apf::ModelEntity* getFaceModelEntity_yz(const int i, const int j, const int k)
    {
      int xmax = m_meshspec->nx;

      ModelEntitySpec me;
      if (i == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[3]);
      else if (i == xmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[1]);
      else
        me = ModelEntitySpec(3, m_current_block);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for faces parallel to the xy plane
    apf::ModelEntity* getFaceModelEntity_xy(const int i, const int j, const int k)
    {
      int zmax = m_meshspec->nz;

      ModelEntitySpec me;
      if (k == 0)
        me = ModelEntitySpec(2, m_block_geometric_faces[4]);
      else if (k == zmax)
        me = ModelEntitySpec(2, m_block_geometric_faces[5]);
      else
        me = ModelEntitySpec(3, m_current_block);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    void createElements()
    {
      for (int block=0; block < m_meshspecs.size(); ++block)
      {
        setBlock(block);
        apf::MeshEntity* verts[8];
        apf::ModelEntity* me = m_mesh->findModelEntity(3, m_current_block);
        for (int i=0; i < m_meshspec->nx; ++i)
          for (int j=0; j < m_meshspec->ny; ++j)
            for (int k=0; k < m_meshspec->nz; ++k)
            {
              int j2 = j + m_yidx_offset;
              verts[0] = m_verts[i][j2][k];
              verts[1] = m_verts[i+1][j2][k];
              verts[2] = m_verts[i+1][j2+1][k];
              verts[3] = m_verts[i][j2+1][k];
              verts[4] = m_verts[i][j2][k+1];
              verts[5] = m_verts[i+1][j2][k+1];
              verts[6] = m_verts[i+1][j2+1][k+1];
              verts[7] = m_verts[i][j2+1][k+1];
              for (int i=0; i < 8; ++i)
              {
                apf::Vector3 pt;
                m_mesh->getPoint(verts[i], 0, pt);
              }
              
              apf::buildElement(m_mesh, me, apf::Mesh::HEX, verts);
            }
      }
    }

    ModelEntitySpec getMESpec(apf::ModelEntity* me)
    {
      return ModelEntitySpec{m_mesh->getModelType(me), m_mesh->getModelTag(me) };
    }

    void createGeometry()
    {
      std::vector<std::vector<int>> verts_to_edges{{0, 3, 8},
                                                   {0, 1, 9},
                                                   {1, 2, 10},
                                                   {2, 3, 11},
                                                   {4, 7, 8},
                                                   {4, 5, 9},
                                                   {5, 6, 10},
                                                   {6, 7, 11}
                                                   };
      std::vector<std::vector<int>> edges_to_faces{{0, 4},
                                                   {1, 4},
                                                   {2, 4},
                                                   {3, 4},
                                                   {0, 5},
                                                   {1, 5},
                                                   {2, 5},
                                                   {3, 5},
                                                   {0, 3},
                                                   {0, 1},
                                                   {1, 2},
                                                   {2, 3}};
      // each element has all the faces, so no need to have a vector for it
      for (unsigned int i=0; i < m_meshspecs.size(); ++i)
      {
        setBlock(i);
        for (int dim=0; dim < 3; ++dim)
          for (unsigned int j=0; j < m_block_geometric_entities[dim].get().size(); ++j)
          {
            int tag = m_block_geometric_entities[dim].get()[j];
            if (!m_gmi_topo->hasEntityByTag(dim, tag))
            {
              m_gmi_topo->createEntity(dim, tag);
            }
          }

        m_gmi_topo->createEntity(3, i);
      }

      for (unsigned int i=0; i < m_meshspecs.size(); ++i)
      {
        setBlock(i);
        for (int dim=0; dim < 3; ++dim)
          for (int j=0; j < m_block_geometric_entities[dim].get().size(); ++j)
          {
            int tag = m_block_geometric_entities[dim].get()[j];
            auto& entity = m_gmi_topo->getEntityByTag(dim, tag);
            if (dim == 2)
            {
              auto& element = m_gmi_topo->getEntityByTag(3, i);
              entity.addUpward(element);
              element.addDownward(entity);
            } else
            {
              auto& upward_block_indices = dim == 0 ? verts_to_edges[j] : edges_to_faces[j];
              createConnectivity(entity, upward_block_indices);
            }
          }
      }

      mesh_gmi::verify(*m_gmi_topo);
    
                                                         
    //TODO: make GMIEntity non-copyable
    }

    bool hasUpwardAdjacency(const mesh_gmi::GMIEntity& entity, const mesh_gmi::GMIEntity& entity_up)
    {
      for (int i=0; i < entity.countUpward(); ++i)
        if (entity.getUpwardIndex(i) == entity_up.getIndex())
          return true;
      
      return false;
    }

    void createConnectivity(mesh_gmi::GMIEntity& entity, const std::vector<int>& upward_block_indices)
    {
      int dim = entity.getDim();
      for (unsigned int i=0; i < upward_block_indices.size(); ++i)
      {
        int tag = m_block_geometric_entities[dim+1].get()[upward_block_indices[i]];
        auto& up_entity = m_gmi_topo->getEntityByTag(dim+1, tag);
        if (!hasUpwardAdjacency(entity, up_entity))
        {
          entity.addUpward(up_entity);
          up_entity.addDownward(entity);
        }
      }
    }


    void setBlock(int idx)
    {
      assert(idx >= 0 && unsigned(idx) < m_meshspecs.size());
      m_current_block = idx;

      int idx_offset = 0;
      Real coord_offset = m_meshspecs[0].ymin;
      for (int i=0; i < idx; ++i)
      {
        idx_offset += m_meshspecs[i].ny;
        coord_offset += m_meshspecs[i].ymax - m_meshspecs[i].ymin;
      }

      m_coords_min[0] = m_meshspecs[0].xmin;
      m_coords_min[1] = coord_offset;
      m_coords_min[2] = m_meshspecs[0].ymin;
      m_yidx_offset   = idx_offset;
      m_meshspec      = &(m_meshspecs[idx]);

      setGeometricVerts(idx);
      setGeometricEdges(idx);
      setGeometricFaces(idx);
    }

    void setGeometricVerts(int block)
    {
      if (block == 0)
      {
        for (int i=0; i < 8; ++i)
          m_block_geometric_verts[i] = i;
      } else if (block == 1)
      {
        std::array<int, 8> new_vals = {3, 2, 8, 9, 7, 6, 10, 11};
        for (int i=0; i < 8; ++i)
          m_block_geometric_verts[i] = new_vals[i];     
      } else
      {
        m_block_geometric_verts[0] = m_block_geometric_verts[3];
        m_block_geometric_verts[3] += 4;

        m_block_geometric_verts[1] = m_block_geometric_verts[2];
        m_block_geometric_verts[2] += 4;

        m_block_geometric_verts[4] = m_block_geometric_verts[7];
        m_block_geometric_verts[7] += 4;

        m_block_geometric_verts[5] = m_block_geometric_verts[6];
        m_block_geometric_verts[6] += 4;
      }
    }

    void setGeometricEdges(int block)
    {
      if (block == 0)
      {
        for (int i=0; i < 12; ++i)
          m_block_geometric_edges[i] = i;
      } else if (block == 1)
      {
        std::array<int, 12> new_vals = {2, 12, 13, 14, 6, 15, 16, 17, 11, 10, 18, 19};
        for (int i=0; i < 12; ++i)
          m_block_geometric_edges[i] = new_vals[i];     
      } else
      {
        std::vector<std::pair<int, int>> shuffles = { {0, 2}, {9, 10}, {8, 11}, {4, 6} };
        int offset = 8;
        for (auto& v : m_block_geometric_edges)
          v += offset;

        for (auto& p : shuffles)
          m_block_geometric_edges[p.first] = m_block_geometric_edges[p.second] - offset;
      }
    }

    void setGeometricFaces(int block)
    {
      if (block == 0)
      {
        for (int i=0; i < 6; ++i)
          m_block_geometric_faces[i] = i;
      } else if (block == 1)
      {
        std::array<int, 6> new_vals = {2, 6, 7, 8, 9, 10};
        for (int i=0; i < 6; ++i)
          m_block_geometric_faces[i] = new_vals[i];
      } else
      {
        m_block_geometric_faces[0] = m_block_geometric_faces[2];
        for (int i=1; i < 6; ++i)
          m_block_geometric_faces[i] += 5;
      }
    }

    void checkInput()
    {
      assertAlways(m_meshspecs.size() > 0, "must have at least 1 MeshSpec for multi-block mesh");
      auto& meshspec0 = m_meshspecs[0];
      for (unsigned int i=1; i < m_meshspecs.size(); ++i)
      {
        auto& meshspec = m_meshspecs[i];
        assertAlways(meshspec0.nx == meshspec.nx, "Number of elements in x direction must be same for all blocks");
        assertAlways(meshspec0.nz == meshspec.nz, "Number of elements in z direction must be same for all blocks");
        assertAlways(std::abs(meshspec0.xmin - meshspec.xmin) < 1e-13, "x dimension must be same for all blocks");
        assertAlways(std::abs(meshspec0.xmax - meshspec.xmax) < 1e-13, "x dimension must be same for all blocks");
        assertAlways(std::abs(meshspec0.zmin - meshspec.zmin) < 1e-13, "z dimension must be same for all blocks");
        assertAlways(std::abs(meshspec0.zmax - meshspec.zmax) < 1e-13, "z dimension must be same for all blocks");

        auto& meshspec_prev = m_meshspecs[i-1];
        assertAlways(std::abs(meshspec_prev.ymax - meshspec.ymin) < 1e-13, "Blocks must be contiguous in y dimension");
      }
    }

    T m_func; // mapping function
    std::vector<MeshSpec> m_meshspecs;
    MeshSpec* m_meshspec;
    ArrayType<apf::MeshEntity*, 3> m_verts;
    std::shared_ptr<mesh_gmi::GMITopo> m_gmi_topo;
    gmi_model* m_gmi_model;
    apf::Mesh2* m_mesh;
    int m_current_block = 0;
    int m_yidx_offset = 0;
    std::array<Real, 3> m_coords_min;
    std::vector<int>  m_block_geometric_verts;
    std::vector<int>  m_block_geometric_edges;
    std::vector<int>  m_block_geometric_faces;
    std::vector<std::reference_wrapper<std::vector<int>>> m_block_geometric_entities;


};

template <typename T>
MeshGeneratorMulti<typename std::decay<T>::type>
make_mesh_generator(const std::vector<MeshSpec>& meshspecs, T&& func=&identity)
{
  return {meshspecs, std::forward<T>(func)};
}

//TODO: do instantiation in source file
extern template class MeshGeneratorMulti<IdentityType>;

} // namespace

#endif  // header guard
