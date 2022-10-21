#ifndef MESH_GEOMETRY_MULTI_H
#define MESH_GEOMETRY_MULTI_H

#include "ProjectDefs.h"
#include "error_handling.h"
#include "mesh/gmiDataStructure.h"
#include "mesh/mesh_input.h"
#include <apfVector.h>
#include <gmi_null.h>
#include <iostream>
#include <vector>
#include <memory>

#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"

namespace apf {
inline std::ostream& operator<<(std::ostream& os, const Vector3& vec)
{
  os << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
  return os;
}
}

namespace Mesh {


// assigns new geometric entity Ids
class GeometricEntityIDGenerator
{
  public:
    int getNewId(int dim)
    {
      int id = m_id[dim]++;
      return id;
    }

  private:
    std::array<int, 4> m_id = {0, 0, 0, 0};
};

using IdGeneratorPtr = std::shared_ptr<GeometricEntityIDGenerator>;

struct LocalGeometricEntity
{
  explicit LocalGeometricEntity(int dim=0, int id=0) :
    dim(dim),
    id(id)
  {}

  int dim;
  int id;
};

// describes the geometric entity that lives on an 
// adjacent block
struct RelativeGeometricEntity
{
  int block_id;
  int local_geometric_id;
};

struct Indices
{
  int i;
  int j;
  int k;
};


class BlockGeometry
{
  public:
    using BlockVector = std::vector<std::shared_ptr<BlockGeometry>>;
    BlockGeometry(BlockVector& surrounding_blocks, IdGeneratorPtr id_gen, std::shared_ptr<mesh_gmi::GMITopo> gmi_model) :
      m_surrounding_blocks(surrounding_blocks),
      m_gmi_model(gmi_model)
    {
      std::array<int, 4> sizes = {8, 12, 6, 1};
      for (int i=0; i < 4; ++i)
      {
        m_geometric_ids[i].resize(sizes[i]);
        m_found_sources[i].resize(sizes[i]);
        m_geometric_sources[i].resize(sizes[i]);
      }

      setSources();
      setGeometricIds(id_gen);
      createGmiEntities();
      setGmiAdjacency();
    }

    int getGeometricEntity(LocalGeometricEntity ge) const
    {
      return m_geometric_ids[ge.dim][ge.id];
    }

    bool hasGeometricSource(LocalGeometricEntity ge) const
    {
      return m_found_sources[ge.dim][ge.id].block_id != -1;
    }

    RelativeGeometricEntity getGeometricSource(LocalGeometricEntity ge) const
    {
      return m_found_sources[ge.dim][ge.id];
    }

  private:

    void setGeometricIds(IdGeneratorPtr id_gen)
    {
      for (int dim=0; dim < 4; ++dim)
        for (int idx=0; idx < m_geometric_ids[dim].size(); ++idx)
          m_geometric_ids[dim][idx] = getOrCreateGeometricId(dim, idx, id_gen);
    }

    int getOrCreateGeometricId(int dim, int idx, IdGeneratorPtr id_gen)
    {
      auto& sources = m_geometric_sources[dim][idx];
      RelativeGeometricEntity found_source = {-1, -1};
      int entity_id = -1;
      for (auto& source : sources)
      {

        auto block = m_surrounding_blocks[source.block_id];
        if (block)
        {
          entity_id = block->getGeometricEntity(LocalGeometricEntity(dim, source.local_geometric_id));
          found_source = source;
        }
      }

      if (entity_id == -1)
        entity_id = id_gen->getNewId(dim);

      m_geometric_ids[dim][idx] = entity_id;
      m_found_sources[dim][idx] = found_source;

      return entity_id;
    }

    void setSources()
    {
      setVertSources();
      setEdgeSources();
      setFaceSources();
      setDomainSources();
    }

    void createGmiEntities()
    {
      createGmiEntities(0, 8);
      createGmiEntities(1, 12);
      createGmiEntities(2, 6);
      createGmiEntities(3, 1);
    }

    void createGmiEntities(int dim, int num_entities)
    {
      for (int i=0; i < num_entities; ++i)
      {
        LocalGeometricEntity local_entity(dim, i);
        if (!hasGeometricSource(local_entity))
        {
          int id = getGeometricEntity(local_entity);
          m_gmi_model->createEntity(dim, id);
        }
      }    
    }

    void setGmiAdjacency()
    {
      std::vector<std::pair<int, int>> edge_verts = { {0, 1}, {1, 2}, {2, 3}, {3, 0},
                                                      {4, 5}, {5, 6}, {6, 7}, {7, 4},
                                                      {0, 4}, {1, 5}, {2, 6}, {3, 7}};
      std::vector<std::array<int, 4>> face_edges = { {0, 1, 2, 3},   {0, 9, 4, 8},  {1, 10, 5, 9},
                                                     {2, 11, 6, 10}, {3, 11, 7, 8}, {4, 5, 6, 7}};

      for (int edge_idx=0; edge_idx < 12; ++edge_idx)
      {
        LocalGeometricEntity edge(1, edge_idx), vert1(0, edge_verts[edge_idx].first),
                             vert2(0, edge_verts[edge_idx].second);
        mesh_gmi::GMIEntity& edge_gmi  = m_gmi_model->getEntityByTag(1, getGeometricEntity(edge));
        mesh_gmi::GMIEntity& vert1_gmi = m_gmi_model->getEntityByTag(0, getGeometricEntity(vert1));
        mesh_gmi::GMIEntity& vert2_gmi = m_gmi_model->getEntityByTag(0, getGeometricEntity(vert2));
        if (!hasGeometricSource(edge) || !hasGeometricSource(vert1))
        {
          edge_gmi.addDownward(vert1_gmi);
          vert1_gmi.addUpward(edge_gmi);
        }

        if (!hasGeometricSource(edge) ||!hasGeometricSource(vert2))
        {
          edge_gmi.addDownward(vert2_gmi);
          vert2_gmi.addUpward(edge_gmi);          
        }
      }

      for (int face_idx=0; face_idx < 6; ++face_idx)
      {
        for (int edge_idx : face_edges[face_idx])
        {
          LocalGeometricEntity face(2, face_idx), edge(1, edge_idx);
          mesh_gmi::GMIEntity& face_gmi  = m_gmi_model->getEntityByTag(2, getGeometricEntity(face));
          mesh_gmi::GMIEntity& edge_gmi  = m_gmi_model->getEntityByTag(1, getGeometricEntity(edge));
          if (!hasGeometricSource(face) || !hasGeometricSource(edge))
          {
            face_gmi.addDownward(edge_gmi);
            edge_gmi.addUpward(face_gmi);
          }
        }
      }

      for (int face_idx=0; face_idx < 6; ++face_idx)
      {
        LocalGeometricEntity domain(3, 0), face(2, face_idx);
        mesh_gmi::GMIEntity& domain_gmi  = m_gmi_model->getEntityByTag(3, getGeometricEntity(domain));
        mesh_gmi::GMIEntity& face_gmi  = m_gmi_model->getEntityByTag(2, getGeometricEntity(face));
        domain_gmi.addDownward(face_gmi);
        face_gmi.addUpward(domain_gmi);
      }
    }

    void setVertSources()
    {
      auto& vert_sources = m_geometric_sources[0];
      vert_sources[0] = { {0, 6}, {1, 7}, {7, 5}, {8, 2},  {9, 3},  {15, 1}, {24, 4} };
      vert_sources[1] = { {1, 6}, {2, 7}, {3, 4}, {9, 2},  {10, 3}, {11, 0}, {24, 5} };
      vert_sources[2] = { {3, 7}, {4, 4}, {5, 5}, {11, 3}, {12, 0}, {13, 1}, {24, 6} };
      vert_sources[3] = { {5, 4}, {6, 5}, {7, 6}, {13, 0}, {14, 1}, {15, 2}, {24, 7} };
      for (int i=0; i < 4; ++i)
      {
        // compute the 6 surrounding entities
        for (int j=0; j < 6; ++j)
        {
          auto& start_entity = vert_sources[i][j];
          RelativeGeometricEntity new_entity = {start_entity.block_id + 8, start_entity.local_geometric_id};
          vert_sources[i+4].push_back(new_entity);
        }

        // compute the entity above
        auto& start_entity = vert_sources[0][6];
        RelativeGeometricEntity new_entity = {start_entity.block_id + 1, start_entity.local_geometric_id - 4};
        vert_sources[i+4].push_back(new_entity);
      }
    }

    void setEdgeSources()
    {
      auto& edge_sources = m_geometric_sources[1];
      edge_sources[0] = { {1, 6}, {9, 2},  {24, 4} };
      edge_sources[1] = { {3, 7}, {11, 3}, {24, 5} };
      edge_sources[2] = { {5, 4}, {13, 0}, {24, 6} };
      edge_sources[3] = { {7, 5}, {15, 1}, {24, 7} };
      std::vector<RelativeGeometricEntity> offsets = { {8, 0}, {8, 0}, {1, -4} };
      for (int i=0; i < 4; ++i)
      {
        // compute the 6 surrounding entities
        for (int j=0; j < 3; ++j)
        {
          auto& start_entity = edge_sources[i][j];
          auto& offset       = offsets[j];
          RelativeGeometricEntity new_entity = {start_entity.block_id + offset.block_id, 
                                                start_entity.local_geometric_id + offset.local_geometric_id};
          edge_sources[i+4].push_back(new_entity);
        }
      }

      edge_sources[8]  = { {8, 10},  {9, 11},  {15, 9}  };
      edge_sources[9]  = { {9, 10},  {10, 11}, {11, 8}  };
      edge_sources[10] = { {11, 11}, {12, 8},  {13, 9}  };
      edge_sources[11] = { {13, 8},  {14, 9},  {15, 10} };
    }

    void setFaceSources()
    {
      auto& face_sources = m_geometric_sources[2];
      face_sources[0] = { {24, 5} };
      face_sources[1] = { {9,  3} };
      face_sources[2] = { {11, 4} };
      face_sources[3] = { {13, 1} };
      face_sources[4] = { {15, 2} };
      face_sources[5] = { {25, 0} };
    }

    void setDomainSources()
    {
      auto& domain_sources = m_geometric_sources[3];
      domain_sources.resize(0);
    }

    std::array< std::vector<int>, 4> m_geometric_ids;
    std::array< std::vector<RelativeGeometricEntity>, 4> m_found_sources;
    BlockVector m_surrounding_blocks;

    //std::vector<int>  m_block_geometric_verts;
    //std::vector<int>  m_block_geometric_edges;
    //std::vector<int>  m_block_geometric_faces;

    using GeometricEntitySources = std::vector< std::vector<RelativeGeometricEntity> >;
    // index as [dim][entity_idx][source_idx]
    std::array<GeometricEntitySources, 4> m_geometric_sources;
    //std::vector< std::vector<RelativeGeometricEntity> > m_vert_sources;
    //std::vector< std::vector<RelativeGeometricEntity> > m_edge_sources;
    //std::vector< std::vector<RelativeGeometricEntity> > m_face_sources;
    std::shared_ptr<mesh_gmi::GMITopo> m_gmi_model;
};


class MeshBlock
{
  public:
    MeshBlock(apf::Mesh2* mesh, MeshSpec meshspec, std::shared_ptr<BlockGeometry> geometry_block, 
              std::vector<std::shared_ptr<MeshBlock>> surrounding_blocks) :
      m_mesh(mesh),
      m_meshspec(meshspec),
      m_coords_min({meshspec.xmin, meshspec.ymin, meshspec.zmin}),
      m_block_geometry(geometry_block),
      m_surrounding_blocks(surrounding_blocks),
      m_verts(boost::extents[meshspec.nx+1][meshspec.ny+1][meshspec.nz+1])
    {}


    void generateMeshBlock()
    {
      createVertices();
      createEdges();
      createFaces();
      createElements();
    }


  private:
    void createVertices()
    {

      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {
            auto local_ge = getVertModelEntity(i, j, k);
            Point p = computeVertCoords(i, j, k);

            if (m_block_geometry->hasGeometricSource(local_ge))
            {
              RelativeGeometricEntity source = m_block_geometry->getGeometricSource(local_ge);
              auto block = m_surrounding_blocks[source.block_id];
              LocalGeometricEntity ge_source(local_ge.dim, source.local_geometric_id);
              m_verts[i][j][k] = block->getMeshVert(ge_source, i, j, k);

              apf::Vector3 p2;
              m_mesh->getPoint(m_verts[i][j][k], 0, p2);
              if (std::abs(p.x - p2.x()) > 1e-13 || std::abs(p.y - p2.y()) > 1e-13 || std::abs(p.z - p2.z()) > 1e-13)
              {
                throw std::runtime_error("retrieve incorrect vert");
              }
            } else
            {
              apf::Vector3 p2(p.x, p.y, p.z);
              apf::ModelEntity* me = m_mesh->findModelEntity(local_ge.dim, m_block_geometry->getGeometricEntity(local_ge));

              m_verts[i][j][k] = m_mesh->createVert(me);
              m_mesh->setPoint(m_verts[i][j][k], 0, p2);
            }
          }
    }


    apf::MeshEntity* getMeshVert(LocalGeometricEntity local_ge, int i, int j, int k)
    {
      int i_max = m_meshspec.nx;
      int j_max = m_meshspec.ny;
      int k_max = m_meshspec.nz;
      std::vector<Indices> vert_indices = { {0, 0, 0},     {i_max, 0, 0},     {i_max, j_max, 0},     {0, j_max, 0},
                                            {0, 0, k_max}, {i_max, 0, k_max}, {i_max, j_max, k_max}, {0, j_max, k_max}};

      std::vector<Indices> edge_indices = { {i, 0, 0},     {i_max, j, 0},     {i, j_max, 0},     {0, j, 0},
                                            {i, 0, k_max}, {i_max, j, k_max}, {i, j_max, k_max}, {0, j, k_max},
                                            {0, 0, k},     {i_max, 0, k},     {i_max, j_max, k}, {0, j_max, k} };
                                         
      std::vector<Indices> face_indices = { {i, j, 0}, {i, 0, k}, {i_max, j, k}, {i, j_max, k}, {0, j, k}, {i, j, k_max} };
      Indices idx;
      if (local_ge.dim == 0)
        idx = vert_indices[local_ge.id];
      else if (local_ge.dim == 1)
        idx = edge_indices[local_ge.id];
      else if (local_ge.dim == 2)
        idx = face_indices[local_ge.id];
      else
        throw std::runtime_error("invalid geometric entity dimension");

      return m_verts[idx.i][idx.j][idx.k];
    }


    LocalGeometricEntity getVertModelEntity(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;
      int ymax = m_meshspec.ny;
      int zmax = m_meshspec.nz;

      LocalGeometricEntity me;
      // findModelEntity(dim, tag)
      // z = zmin vertices
      if (i == 0 && j == 0 && k == 0)
        me = LocalGeometricEntity(0, 0);
      else if (i == xmax && j == 0 && k == 0)
        me = LocalGeometricEntity(0, 1);
      else if (i == xmax && j == ymax && k == 0)
        me = LocalGeometricEntity(0, 2);
      else if (i == 0 && j == ymax && k == 0)
        me = LocalGeometricEntity(0, 3);
      // z = zmax vertices
      else if (i == 0 && j == 0 && k == zmax)
        me = LocalGeometricEntity(0, 4);
      else if (i == xmax && j == 0 && k == zmax)
        me = LocalGeometricEntity(0, 5);
      else if (i == xmax && j == ymax && k == zmax)
        me = LocalGeometricEntity(0, 6);
      else if (i == 0 && j == ymax && k == zmax)
        me = LocalGeometricEntity(0, 7);
      // edges on z = zmin quad
      else if (j == 0 && k == 0)
        me = LocalGeometricEntity(1, 0);
      else if (i == xmax && k == 0)
        me = LocalGeometricEntity(1, 1);
      else if (j == ymax && k == 0)
        me = LocalGeometricEntity(1, 2);
      else if (i == 0 && k == 0)
        me = LocalGeometricEntity(1, 3);
      // edges on z = zmax quad
      else if (j == 0 && k == zmax)
        me = LocalGeometricEntity(1, 4);
      else if (i == xmax && k == zmax)
        me = LocalGeometricEntity(1, 5);
      else if (j == ymax && k == zmax)
        me = LocalGeometricEntity(1, 6);
      else if (i == 0 && k == zmax)
        me = LocalGeometricEntity(1, 7);
      // edges parallel to z axis
      else if (i == 0 && j == 0)
        me = LocalGeometricEntity(1, 8);
      else if (i == xmax && j == 0)
        me = LocalGeometricEntity(1, 9);
      else if (i == xmax && j == ymax)
        me = LocalGeometricEntity(1, 10);
      else if (i == 0 && j == ymax)
        me = LocalGeometricEntity(1, 11);
      // faces
      else if (k == 0)
        me = LocalGeometricEntity(2, 0);
      else if (j == 0)
        me = LocalGeometricEntity(2, 1);
      else if (i == xmax)
        me = LocalGeometricEntity(2, 2);
      else if (j == ymax)
        me = LocalGeometricEntity(2, 3);
      else if (i == 0)
        me = LocalGeometricEntity(2, 4);
      else if (k == zmax)
        me = LocalGeometricEntity(2, 5);
      // interior
      else
        me = LocalGeometricEntity(3, 0);

      return me;
    }


    Point computeVertCoords(const int i, const int j, const int k)
    {
      double dx = (m_meshspec.xmax - m_meshspec.xmin)/m_meshspec.nx;
      double dy = (m_meshspec.ymax - m_meshspec.ymin)/m_meshspec.ny;
      double dz = (m_meshspec.zmax - m_meshspec.zmin)/m_meshspec.nz;

      Point p{i*dx + m_coords_min[0], j*dy + m_coords_min[1], k*dz + m_coords_min[2]};
      return p; //TODO: mapping function
      //return m_func(p);
    }


  void createEdges()
    {
      std::cout << "\nMaking x direction edges" << std::endl;
      apf::MeshEntity* verts[2];
      for (int i=0; i < m_meshspec.nx; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {
            // x direction
            LocalGeometricEntity local_ge = getEdgeModelEntity_x(i, j, k);
            if (m_block_geometry->hasGeometricSource(local_ge))
              continue;

            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i+1][j][k];
            apf::ModelEntity* me = m_mesh->findModelEntity(local_ge.dim, m_block_geometry->getGeometricEntity(local_ge));
            m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
          }

      std::cout << "\nMaking y direction edges" << std::endl;
      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {

            // y direction
            LocalGeometricEntity local_ge = getEdgeModelEntity_y(i, j, k);

            if (m_block_geometry->hasGeometricSource(local_ge))
              continue;

            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i][j+1][k];
            apf::ModelEntity* me = m_mesh->findModelEntity(local_ge.dim, m_block_geometry->getGeometricEntity(local_ge));
            m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
          }

      std::cout << "\nMaking z direction edges" << std::endl;
      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz; ++k)
          {

            // z direction
            LocalGeometricEntity local_ge = getEdgeModelEntity_z(i, j, k);

            if (m_block_geometry->hasGeometricSource(local_ge))
              continue;

            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i][j][k+1];
            apf::ModelEntity* me = m_mesh->findModelEntity(local_ge.dim, m_block_geometry->getGeometricEntity(local_ge));
            m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
          }

    }

    // gets model entity for edges parallel to the x axis
    LocalGeometricEntity getEdgeModelEntity_x(const int i, const int j, const int k)
    {
      //int xmax = m_meshspec.nx + 1;
      int ymax = m_meshspec.ny;
      int zmax = m_meshspec.nz;

      LocalGeometricEntity me;
      // geometric edges
      if (j == 0 && k == 0)
        me = LocalGeometricEntity(1, 0);
      else if (j == ymax && k == 0)
        me = LocalGeometricEntity(1, 2);
      else if (j == 0 && k == zmax)
        me = LocalGeometricEntity(1, 4);
      else if (j == ymax && k == zmax)
        me = LocalGeometricEntity(1, 6);
      // faces
      else if (j == 0)
        me = LocalGeometricEntity(2, 1);
      else if (k == 0)
        me = LocalGeometricEntity(2, 0);
      else if (j == ymax)
        me = LocalGeometricEntity(2, 3);
      else if (k == zmax)
        me = LocalGeometricEntity(2, 5);
      // interior
      else
        me = LocalGeometricEntity(3, 0);


      return me;
    }

    // gets model entity for edges parallel to the y axis
    LocalGeometricEntity getEdgeModelEntity_y(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;
      //int ymax = m_meshspec.ny;
      int zmax = m_meshspec.nz;

      LocalGeometricEntity me;
      // geometric edges
      if (i == 0 && k == 0)
        me = LocalGeometricEntity(1, 3);
      else if (i == xmax && k == 0)
        me = LocalGeometricEntity(1, 1);
      else if (i == 0 && k == zmax)
        me = LocalGeometricEntity(1, 7);
      else if (i == xmax && k == zmax)
        me = LocalGeometricEntity(1, 5);
      // faces
      else if (i == 0)
        me = LocalGeometricEntity(2, 4);
      else if (k == 0)
        me = LocalGeometricEntity(2, 0);
      else if (i == xmax)
        me = LocalGeometricEntity(2, 2);
      else if (k == zmax)
        me = LocalGeometricEntity(2, 5);
      // interior
      else
        me = LocalGeometricEntity(3, 0);
 
      return me;
    }

    // gets model entity for edges parallel to the z axis
    LocalGeometricEntity getEdgeModelEntity_z(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;
      int ymax = m_meshspec.ny;
      //int zmax = m_meshspec.nz;

      LocalGeometricEntity me;
      // geometric edges
      if (j == 0 && i == 0)
        me = LocalGeometricEntity(1, 8);
      else if (j == ymax && i == 0)
        me = LocalGeometricEntity(1, 11);
      else if (j == 0 && i == xmax)
        me = LocalGeometricEntity(1, 9);
      else if (j == ymax && i == xmax)
        me = LocalGeometricEntity(1, 10);
      // faces
      else if (j == 0)
        me = LocalGeometricEntity(2, 1);
      else if (i == 0)
        me = LocalGeometricEntity(2, 4);
      else if (j == ymax)
        me = LocalGeometricEntity(2, 3);
      else if (i == xmax)
        me = LocalGeometricEntity(2, 2);
      // interior
      else
        me = LocalGeometricEntity(3, 0);

      return me;
    }

    apf::MeshEntity* getCommonEdge(apf::MeshEntity* vert1, apf::MeshEntity* vert2)
    {
        for (int j=0; j < m_mesh->countUpward(vert1); ++j)
          for (int k=0; k < m_mesh->countUpward(vert2); ++k)
            if (m_mesh->getUpward(vert1, j) == m_mesh->getUpward(vert2, k))
              return m_mesh->getUpward(vert1, j);

        throw std::runtime_error("unable to find common edge");
    }

    void printEdges(apf::MeshEntity* verts[4])
    {
      for (int i=0; i < 4; ++i)
      {
        int i2 = (i + 1) % 4;
        apf::MeshEntity* edge = getCommonEdge(verts[i], verts[i2]);
        apf::ModelEntity* me = m_mesh->toModel(edge);
        std::cout << "edge " << i << " = " << edge << ", classified on geometric entity " << m_mesh->getModelType(me) << ", " << m_mesh->getModelTag(me) << std::endl;
      }
    }

    void printVerts(apf::MeshEntity* verts[4])
    {
      std::cout << "verts = " << std::endl;
      for (int i=0; i < 4; ++i)
        std::cout << verts[i] << ", ";
      std::cout << std::endl;
    }

    void createFaces()
    {
      apf::MeshEntity* verts[4];
      std::cout << "\ncreating xz faces" << std::endl;
      for (int i=0; i < m_meshspec.nx; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz; ++k)
          {
            // xz plane
            LocalGeometricEntity local_ge = getFaceModelEntity_xz(i, j, k);
            if (m_block_geometry->hasGeometricSource(local_ge))
              continue;   

            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i+1][j][k];
            verts[2] = m_verts[i+1][j][k+1];
            verts[3] = m_verts[i][j][k+1];
            //printEdges(verts);
            apf::ModelEntity* me = m_mesh->findModelEntity(local_ge.dim, m_block_geometry->getGeometricEntity(local_ge));
            apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
          }

      std::cout << "\ncreating yz faces" << std::endl;
      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny; ++j)
          for (int k=0; k < m_meshspec.nz; ++k)
          {
            // yz plane
            LocalGeometricEntity local_ge = getFaceModelEntity_yz(i, j, k);
            //if (m_block_geometry->hasGeometricSource(local_ge))
            //  continue;   

            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i][j][k+1];
            verts[2] = m_verts[i][j+1][k+1];
            verts[3] = m_verts[i][j+1][k];
            //printEdges(verts);
            //printVerts(verts);

            if (m_block_geometry->hasGeometricSource(local_ge))
              continue;   

            apf::ModelEntity* me = m_mesh->findModelEntity(local_ge.dim, m_block_geometry->getGeometricEntity(local_ge));
            apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
          }

      std::cout << "\ncreating xy faces" << std::endl;
      for (int i=0; i < m_meshspec.nx; ++i)
        for (int j=0; j < m_meshspec.ny; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {

            // xy plane
            LocalGeometricEntity local_ge = getFaceModelEntity_xy(i, j, k);
            if (m_block_geometry->hasGeometricSource(local_ge))
              continue;  

            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i+1][j][k];
            verts[2] = m_verts[i+1][j+1][k];
            verts[3] = m_verts[i][j+1][k];
            //printEdges(verts);
            apf::ModelEntity* me = m_mesh->findModelEntity(local_ge.dim, m_block_geometry->getGeometricEntity(local_ge));            
            apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
          }
    }

    // gets model entity for faces parallel to the xz plane
    LocalGeometricEntity getFaceModelEntity_xz(const int i, const int j, const int k)
    {
      int ymax = m_meshspec.ny;

      LocalGeometricEntity me;
      if (j == 0)
        me = LocalGeometricEntity(2, 1);
      else if (j == ymax)
        me = LocalGeometricEntity(2, 3);
      else
        me = LocalGeometricEntity(3, 0);

      return me;
    }

    // gets model entity for faces parallel to the yz plane
    LocalGeometricEntity getFaceModelEntity_yz(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;

      LocalGeometricEntity me;
      if (i == 0)
        me = LocalGeometricEntity(2, 4);
      else if (i == xmax)
        me = LocalGeometricEntity(2, 2);
      else
        me = LocalGeometricEntity(3, 0);

      return me;
    }

    // gets model entity for faces parallel to the xy plane
    LocalGeometricEntity getFaceModelEntity_xy(const int i, const int j, const int k)
    {
      int zmax = m_meshspec.nz;

      LocalGeometricEntity me;
      if (k == 0)
        me = LocalGeometricEntity(2, 0);
      else if (k == zmax)
        me = LocalGeometricEntity(2, 5);
      else
        me = LocalGeometricEntity(3, 0);

      return me;
    }

    void createElements()
    {
      std::cout << "\nCreating elements" << std::endl;
      apf::MeshEntity* verts[8];
      apf::ModelEntity* me = m_mesh->findModelEntity(3, m_block_geometry->getGeometricEntity(LocalGeometricEntity(3, 0)));
      for (int i=0; i < m_meshspec.nx; ++i)
        for (int j=0; j < m_meshspec.ny; ++j)
          for (int k=0; k < m_meshspec.nz; ++k)
          {
            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i+1][j][k];
            verts[2] = m_verts[i+1][j+1][k];
            verts[3] = m_verts[i][j+1][k];
            verts[4] = m_verts[i][j][k+1];
            verts[5] = m_verts[i+1][j][k+1];
            verts[6] = m_verts[i+1][j+1][k+1];
            verts[7] = m_verts[i][j+1][k+1];
            
            apf::buildElement(m_mesh, me, apf::Mesh::HEX, verts);
          }
    }


    apf::Mesh2* m_mesh;
    MeshSpec m_meshspec;
    std::array<Real, 3> m_coords_min;
    std::shared_ptr<BlockGeometry> m_block_geometry;
    std::vector<std::shared_ptr<MeshBlock>> m_surrounding_blocks;
    ArrayType<apf::MeshEntity*, 3> m_verts;
};


struct MultiBlockMeshSpec
{
  MeshSpec middle_block;  // Note: this is required, even if the block is not actually
                          //       created, because some dimensions of the other blocks
                          //       are inferred from it
  bool create_middle_block;

  // Note: the element counts and thicknesses are defined outward
  //       from the middle block
  std::vector<int> numel_plusx;
  std::vector<Real> thickness_plusx;

  std::vector<int> numel_minusx;
  std::vector<Real> thickness_minusx;  

  std::vector<int> numel_plusy;
  std::vector<Real> thickness_plusy;
  std::vector<int> numel_minusy;
  std::vector<Real> thickness_minusy;  

  std::vector<int> numel_plusz;
  std::vector<Real> thickness_plusz;

  std::vector<int> numel_minusz;
  std::vector<Real> thickness_minusz;
};

void validateMultiBlockMeshSpec(const MultiBlockMeshSpec& spec)
{
  if (spec.numel_plusx.size() != spec.thickness_plusx.size())
    throw std::runtime_error("plusx number of layers is inconsistent");

  if (spec.numel_plusy.size() != spec.thickness_plusy.size())
      throw std::runtime_error("plusy number of layers is inconsistent");    

  if (spec.numel_plusz.size() != spec.thickness_plusz.size())
      throw std::runtime_error("plusz number of layers is inconsistent");

  if (spec.numel_minusx.size() != spec.thickness_minusx.size())
    throw std::runtime_error("minusx number of layers is inconsistent");

  if (spec.numel_minusy.size() != spec.thickness_minusy.size())
      throw std::runtime_error("minusy number of layers is inconsistent");    

  if (spec.numel_minusz.size() != spec.thickness_minusz.size())
      throw std::runtime_error("minusz number of layers is inconsistent");    

}


//TODO: this can go in source file
template <typename T>
std::vector<T> collectValues(const std::vector<T>& values_minus, T value_middle, 
                             const std::vector<T>& values_plus)
{
  std::vector<T> values;
  for (auto it = values_minus.rbegin(); it != values_minus.rend(); ++it)
    values.push_back(*it);

  values.push_back(value_middle);

  for (auto thickness : values_plus)
    values.push_back(thickness);

  return values;
}

class MeshGeneratorMultiBlock
{
  public:
    MeshGeneratorMultiBlock(MultiBlockMeshSpec meshspec) :
      m_meshspec(meshspec),
      m_geometric_id_gen(std::make_shared<GeometricEntityIDGenerator>())
    {
      validateMultiBlockMeshSpec(meshspec);
      reformatBlockData();
      computeBlockMeshSpecs();
    }


    apf::Mesh2* generate()
    {
      initializeGeometry();
      m_mesh = apf::makeEmptyMdsMesh(m_gmi_model, 3, false);
      createMeshBlocks();

      m_mesh->acceptChanges();
      m_mesh->verify();

      apf::writeASCIIVtkFiles("mesh_created", m_mesh);
     
      return m_mesh;
    }

    std::shared_ptr<mesh_gmi::GMITopo> getGmiTopo() { return m_gmi_topo; }

  private:

    void reformatBlockData()
    {
      m_thicknesses[0] = collectValues(m_meshspec.thickness_minusx, 
                                       m_meshspec.middle_block.xmax - m_meshspec.middle_block.xmin,
                                       m_meshspec.thickness_plusx);
      m_thicknesses[1] = collectValues(m_meshspec.thickness_minusy, 
                                       m_meshspec.middle_block.ymax - m_meshspec.middle_block.ymin,
                                       m_meshspec.thickness_plusy);

      m_thicknesses[2] = collectValues(m_meshspec.thickness_minusz, 
                                       m_meshspec.middle_block.zmax - m_meshspec.middle_block.zmin,
                                       m_meshspec.thickness_plusz);

      m_numels[0] = collectValues(m_meshspec.numel_minusx, m_meshspec.middle_block.nx, m_meshspec.numel_plusx);
      m_numels[1] = collectValues(m_meshspec.numel_minusy, m_meshspec.middle_block.ny, m_meshspec.numel_plusy);
      m_numels[2] = collectValues(m_meshspec.numel_minusz, m_meshspec.middle_block.nz, m_meshspec.numel_plusz);
    }

    void computeBlockMeshSpecs()
    {
      int nblocks_x = m_numels[0].size();
      int nblocks_y = m_numels[1].size();
      int nblocks_z = m_numels[2].size();

      m_lower_corner_coords = getLowerCornerCoords();
      m_meshspecs.resize(boost::extents[nblocks_x][nblocks_y][nblocks_z]);
      std::array<Real, 3> current_lower_corner_coords = m_lower_corner_coords;
      for (int i=0; i < nblocks_x; ++i)
      {
        current_lower_corner_coords[1] = m_lower_corner_coords[1];
        for (int j=0; j < nblocks_y; ++j)
        {
          current_lower_corner_coords[2] = m_lower_corner_coords[2];
          for (int k=0; k < nblocks_z; ++k)
          {
            MeshSpec spec_ijk;
            spec_ijk.xmin = current_lower_corner_coords[0];
            spec_ijk.ymin = current_lower_corner_coords[1];
            spec_ijk.zmin = current_lower_corner_coords[2];

            spec_ijk.xmax = spec_ijk.xmin + m_thicknesses[0][i];
            spec_ijk.ymax = spec_ijk.ymin + m_thicknesses[1][j];
            spec_ijk.zmax = spec_ijk.zmin + m_thicknesses[2][k];

            spec_ijk.nx = m_numels[0][i];
            spec_ijk.ny = m_numels[1][j];
            spec_ijk.nz = m_numels[2][k];

            m_meshspecs[i][j][k] = spec_ijk;
            std::cout << "block " << i << ", " << j << ", " << k << " meshspec = " << spec_ijk << std::endl;

            current_lower_corner_coords[2] += m_thicknesses[2][k];
          }

          current_lower_corner_coords[1] += m_thicknesses[1][j];
        }
        current_lower_corner_coords[0] += m_thicknesses[0][i];
      }
    }

    std::array<Real, 3> getLowerCornerCoords()
    {
      std::array<Real, 3> coords = {m_meshspec.middle_block.xmin, m_meshspec.middle_block.ymin,
                                    m_meshspec.middle_block.zmin};

      for (auto thickness : m_meshspec.thickness_minusx)
        coords[0] -= thickness;

      for (auto thickness : m_meshspec.thickness_minusy)
        coords[1] -= thickness;

      for (auto thickness : m_meshspec.thickness_minusz)
        coords[2] -= thickness;

      return coords;
    }

    void initializeGeometry()
    {
      gmi_register_null();
      m_gmi_topo = std::make_shared<mesh_gmi::GMITopo>();
      createGeometryBlocks();
      m_gmi_model = mesh_gmi::createGMITopo(m_gmi_topo);
      //m_gmi_model = gmi_load(".null");
    }

    void createGeometryBlocks()
    {
      m_geometry_blocks.resize(boost::extents[m_numels[0].size()][m_numels[1].size()][m_numels[2].size()]);
      for (int i=0; i < m_numels[0].size(); ++i)
        for (int j=0; j < m_numels[1].size(); ++j)
          for (int k=0; k < m_numels[2].size(); ++k)
          {
            std::cout << "\ncreating geometric block " << i << ", " << j << ", " << k << std::endl;
            auto surrounding_blocks = getSurroundingGeometryBlocks(i, j, k);
            m_geometry_blocks[i][j][k] = std::make_shared<BlockGeometry>(surrounding_blocks, m_geometric_id_gen, m_gmi_topo);
          }
    }


    void createMeshBlocks()
    {
      m_mesh_blocks.resize(boost::extents[m_numels[0].size()][m_numels[1].size()][m_numels[2].size()]);
      for (int i=0; i < m_numels[0].size(); ++i)
        for (int j=0; j < m_numels[1].size(); ++j)
          for (int k=0; k < m_numels[2].size(); ++k)
          {
            std::cout << "\ncreating mesh block " << i << ", " << j << ", " << k << std::endl;

            auto surrounding_blocks = getSurroundingMeshBlocks(i, j, k);
            m_mesh_blocks[i][j][k] = std::make_shared<MeshBlock>(m_mesh, m_meshspecs[i][j][k],
                                       m_geometry_blocks[i][j][k], surrounding_blocks);
            m_mesh_blocks[i][j][k]->generateMeshBlock();
          }
    }

    std::vector<std::shared_ptr<BlockGeometry>> getSurroundingGeometryBlocks(int i, int j, int k)
    {
      std::vector<std::shared_ptr<BlockGeometry>> blocks;
      std::vector<Indices> zplane_indices = { {i-1, j-1, 0}, {i, j-1, 0}, {i+1, j-1, 0},
                                              {i+1, j, 0},   {i+1, j+1, 0},
                                              {i, j+1, 0},   {i-1, j+1, 0},
                                              {i-1, j}};
      std::shared_ptr<BlockGeometry> block;
      for (int zidx=k-1; zidx <= k+1; ++zidx)
      {
        for (auto idx : zplane_indices)
        {
          block = isBlockValid(idx.i, idx.j, zidx) ? m_geometry_blocks[idx.i][idx.j][zidx] : nullptr;
          blocks.push_back(block);
        }
      }

      block = isBlockValid(i, j, k-1) ? m_geometry_blocks[i][j][k-1] : nullptr;
      blocks.push_back(block);

      block = isBlockValid(i, j, k+1) ? m_geometry_blocks[i][j][k+1] : nullptr;
      blocks.push_back(block);

      return blocks;
    }

    std::vector<std::shared_ptr<MeshBlock>> getSurroundingMeshBlocks(int i, int j, int k)
    {
      std::vector<std::shared_ptr<MeshBlock>> blocks;
      std::vector<Indices> zplane_indices = { {i-1, j-1, 0}, {i, j-1, 0}, {i+1, j-1, 0},
                                              {i+1, j, 0},   {i+1, j+1, 0},
                                              {i, j+1, 0},   {i-1, j+1, 0},
                                              {i-1, j}};
      std::shared_ptr<MeshBlock> block;
      for (int zidx=k-1; zidx <= k+1; ++zidx)
      {
        for (auto idx : zplane_indices)
        {
          block = isBlockValid(idx.i, idx.j, zidx) ? m_mesh_blocks[idx.i][idx.j][zidx] : nullptr;
          blocks.push_back(block);
        }
      }

      block = isBlockValid(i, j, k-1) ? m_mesh_blocks[i][j][k-1] : nullptr;
      blocks.push_back(block);

      block = isBlockValid(i, j, k+1) ? m_mesh_blocks[i][j][k+1] : nullptr;
      blocks.push_back(block);
      


      return blocks;
    }

    bool isBlockValid(int i, int j, int k)
    {
      int imiddle = m_meshspec.numel_minusx.size();
      int jmiddle = m_meshspec.numel_minusy.size();
      int kmiddle = m_meshspec.numel_minusz.size();
      if (i == imiddle && j == jmiddle && k == kmiddle)
        return m_meshspec.create_middle_block;
      else
        return i >= 0 && i < m_numels[0].size() &&
               j >= 0 && j < m_numels[1].size() &&
               k >= 0 && k < m_numels[2].size();
    }



    MultiBlockMeshSpec m_meshspec;
    ArrayType<MeshSpec, 3> m_meshspecs;
    std::array<std::vector<Real>, 3> m_thicknesses;
    std::array<std::vector<int>, 3> m_numels;
    std::array<Real, 3> m_lower_corner_coords;

    std::shared_ptr<mesh_gmi::GMITopo> m_gmi_topo;
    gmi_model* m_gmi_model;
    apf::Mesh2* m_mesh;
    std::shared_ptr<GeometricEntityIDGenerator> m_geometric_id_gen;
    ArrayType<std::shared_ptr<BlockGeometry>, 3> m_geometry_blocks;
    ArrayType<std::shared_ptr<MeshBlock>, 3> m_mesh_blocks;
};




}  // namespace

#endif