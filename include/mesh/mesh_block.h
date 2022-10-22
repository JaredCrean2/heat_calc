#ifndef MESH_GENERATOR_MESH_BLOCK_H
#define MESH_GENERATOR_MESH_BLOCK_H

#include "mesh_geometry_multi.h"

namespace Mesh {

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


    void generateMeshBlock();

  private:
    void createVertices();


    apf::MeshEntity* getMeshVert(LocalGeometricEntity local_ge, int i, int j, int k);


    LocalGeometricEntity getVertModelEntity(const int i, const int j, const int k);


    Point computeVertCoords(const int i, const int j, const int k);


    void createEdges();

    // gets model entity for edges parallel to the x axis
    LocalGeometricEntity getEdgeModelEntity_x(const int i, const int j, const int k);

    // gets model entity for edges parallel to the y axis
    LocalGeometricEntity getEdgeModelEntity_y(const int i, const int j, const int k);

    // gets model entity for edges parallel to the z axis
    LocalGeometricEntity getEdgeModelEntity_z(const int i, const int j, const int k);

    apf::MeshEntity* getCommonEdge(apf::MeshEntity* vert1, apf::MeshEntity* vert2);

    void printEdges(apf::MeshEntity* verts[4]);

    void printVerts(apf::MeshEntity* verts[4]);

    void createFaces();

    // gets model entity for faces parallel to the xz plane
    LocalGeometricEntity getFaceModelEntity_xz(const int i, const int j, const int k);

    // gets model entity for faces parallel to the yz plane
    LocalGeometricEntity getFaceModelEntity_yz(const int i, const int j, const int k);

    // gets model entity for faces parallel to the xy plane
    LocalGeometricEntity getFaceModelEntity_xy(const int i, const int j, const int k);

    void createElements();


    apf::Mesh2* m_mesh;
    MeshSpec m_meshspec;
    std::array<Real, 3> m_coords_min;
    std::shared_ptr<BlockGeometry> m_block_geometry;
    std::vector<std::shared_ptr<MeshBlock>> m_surrounding_blocks;
    ArrayType<apf::MeshEntity*, 3> m_verts;
};

} // namespace

#endif