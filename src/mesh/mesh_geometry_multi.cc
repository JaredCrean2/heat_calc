#include "mesh/mesh_geometry_multi.h"

namespace Mesh {


BlockGeometry::BlockGeometry(BlockVector& surrounding_blocks, IdGeneratorPtr id_gen, std::shared_ptr<mesh_gmi::GMITopo> gmi_model) :
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

void BlockGeometry::setGeometricIds(IdGeneratorPtr id_gen)
{
  for (int dim=0; dim < 4; ++dim)
    for (int idx=0; idx < m_geometric_ids[dim].size(); ++idx)
      m_geometric_ids[dim][idx] = getOrCreateGeometricId(dim, idx, id_gen);
}

int BlockGeometry::getOrCreateGeometricId(int dim, int idx, IdGeneratorPtr id_gen)
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

void BlockGeometry::setSources()
{
  setVertSources();
  setEdgeSources();
  setFaceSources();
  setDomainSources();
}

void BlockGeometry::createGmiEntities()
{
  createGmiEntities(0, 8);
  createGmiEntities(1, 12);
  createGmiEntities(2, 6);
  createGmiEntities(3, 1);
}

void BlockGeometry::createGmiEntities(int dim, int num_entities)
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

void BlockGeometry::setGmiAdjacency()
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

void BlockGeometry::setVertSources()
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

void BlockGeometry::setEdgeSources()
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

void BlockGeometry::setFaceSources()
{
  auto& face_sources = m_geometric_sources[2];
  face_sources[0] = { {24, 5} };
  face_sources[1] = { {9,  3} };
  face_sources[2] = { {11, 4} };
  face_sources[3] = { {13, 1} };
  face_sources[4] = { {15, 2} };
  face_sources[5] = { {25, 0} };
}

void BlockGeometry::setDomainSources()
{
  auto& domain_sources = m_geometric_sources[3];
  domain_sources.resize(0);
}


}  // namespace