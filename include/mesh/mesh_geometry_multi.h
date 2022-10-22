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
    BlockGeometry(BlockVector& surrounding_blocks, IdGeneratorPtr id_gen, std::shared_ptr<mesh_gmi::GMITopo> gmi_model);

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

    void setGeometricIds(IdGeneratorPtr id_gen);

    int getOrCreateGeometricId(int dim, int idx, IdGeneratorPtr id_gen);

    void setSources();

    void createGmiEntities();

    void createGmiEntities(int dim, int num_entities);

    void setGmiAdjacency();

    void setVertSources();

    void setEdgeSources();

    void setFaceSources();

    void setDomainSources();

    std::array< std::vector<int>, 4> m_geometric_ids;
    std::array< std::vector<RelativeGeometricEntity>, 4> m_found_sources;
    BlockVector m_surrounding_blocks;

    using GeometricEntitySources = std::vector< std::vector<RelativeGeometricEntity> >;
    std::array<GeometricEntitySources, 4> m_geometric_sources;

    std::shared_ptr<mesh_gmi::GMITopo> m_gmi_model;
};





}  // namespace

#endif