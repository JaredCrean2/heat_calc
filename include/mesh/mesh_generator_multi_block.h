#ifndef MESH_GENERATOR_MULTI_BLOCK_H
#define MESH_GENERATOR_MULTI_BLOCK_H

#include "mesh_block.h"

namespace Mesh {

struct MultiBlockMeshSpec
{
  MeshSpec middle_block;  // Note: this is required, even if the block is not actually
                          //       created, because some dimensions of the other blocks
                          //       are inferred from it

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

  ArrayType<bool, 3> create_blocks;
};

void validateMultiBlockMeshSpec(MultiBlockMeshSpec& spec);


class MeshGeneratorMultiBlock
{
  public:
    MeshGeneratorMultiBlock(MultiBlockMeshSpec meshspec);

    apf::Mesh2* generate();

    std::shared_ptr<mesh_gmi::GMITopo> getGmiTopo();

    int getSurfaceGeometricId(int i, int j, int k, int face);

    int getVolumeGeometricId(int i, int j, int k);

    // indices are in the "middle block" coordinate system: the middle block has
    // indices (0, 0, 0)
    std::shared_ptr<BlockGeometry> getGeometryBlock(int i, int j, int k);

  private:

    void reformatBlockData();

    void computeBlockMeshSpecs();

    std::array<Real, 3> getLowerCornerCoords();

    void initializeGeometry();

    void createGeometryBlocks();

    void createMeshBlocks();

    std::vector<std::shared_ptr<BlockGeometry>> getSurroundingGeometryBlocks(int i, int j, int k);

    std::vector<std::shared_ptr<MeshBlock>> getSurroundingMeshBlocks(int i, int j, int k);

    bool isBlockValid(int i, int j, int k);



    MultiBlockMeshSpec m_meshspec;
    ArrayType<MeshSpec, 3> m_meshspecs;
    std::array<std::vector<Real>, 3> m_thicknesses;
    std::array<std::vector<int>, 3> m_numels;
    std::array<Real, 3> m_lower_corner_coords;
    std::array<int, 3> m_middle_block_indices;

    std::shared_ptr<mesh_gmi::GMITopo> m_gmi_topo;
    gmi_model* m_gmi_model;
    apf::Mesh2* m_mesh;
    std::shared_ptr<GeometricEntityIDGenerator> m_geometric_id_gen;
    ArrayType<std::shared_ptr<BlockGeometry>, 3> m_geometry_blocks;
    ArrayType<std::shared_ptr<MeshBlock>, 3> m_mesh_blocks;
};

}  // namespace

#endif
