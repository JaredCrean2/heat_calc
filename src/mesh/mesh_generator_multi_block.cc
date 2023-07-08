#include "mesh/mesh_generator_multi_block.h"
#include "mesh/apfSplit.h"
#include "utils/mpi_utils.h"
#include <PCU.h>

namespace Mesh {

void validateMultiBlockMeshSpec(MultiBlockMeshSpec& spec)
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

  int nx = 1 + spec.numel_plusx.size() + spec.numel_minusx.size();
  int ny = 1 + spec.numel_plusy.size() + spec.numel_minusy.size();
  int nz = 1 + spec.numel_plusz.size() + spec.numel_minusz.size();
  if (spec.create_blocks.shape()[0] == 0 && spec.create_blocks.shape()[1] == 0 &&
      spec.create_blocks.shape()[2] == 0)
  {
    spec.create_blocks.resize(boost::extents[nx][ny][nz]);
    for (int i=0; i < nx; ++i)
      for (int j=0; j < ny; ++j)
        for (int k=0; k < nz; ++k)
          spec.create_blocks[i][j][k] = true;
  } else
  {
    if (spec.create_blocks.shape()[0] != nx)
      throw std::runtime_error("create_blocks x dimension is incorrect");

    if (spec.create_blocks.shape()[1] != ny)
      throw std::runtime_error("create_blocks y dimension is incorrect");

    if (spec.create_blocks.shape()[2] != nz)
      throw std::runtime_error("create_blocks z dimension is incorrect");
  }
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

MeshGeneratorMultiBlock::MeshGeneratorMultiBlock(MultiBlockMeshSpec meshspec) :
  m_meshspec(meshspec),
  m_geometric_id_gen(std::make_shared<GeometricEntityIDGenerator>())
{
  validateMultiBlockMeshSpec(m_meshspec);
  reformatBlockData();
  computeBlockMeshSpecs();
}


apf::Mesh2* MeshGeneratorMultiBlock::generate()
{
  initializeGeometry();
  m_mesh = nullptr;
  
  MPI_Comm comm_orig = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_SELF);
  if (commRank(comm_orig) == 0)
  {    
    m_mesh = apf::makeEmptyMdsMesh(m_gmi_model, 3, false);
    createMeshBlocks();
    m_mesh->acceptChanges();
    m_mesh->verify();    
  }


  mesh::ApfSplitter splitter(m_mesh, m_gmi_model);
  apf::Mesh2* mesh_parallel = splitter.split(commSize(comm_orig), comm_orig);
  //apf::writeASCIIVtkFiles("mesh_split", m_mesh);

  mesh_parallel->acceptChanges();
  mesh_parallel->verify();

  apf::writeASCIIVtkFiles("mesh_created", mesh_parallel);
  
  return mesh_parallel;
}

std::shared_ptr<mesh_gmi::GMITopo> MeshGeneratorMultiBlock::getGmiTopo() { return m_gmi_topo; }

int MeshGeneratorMultiBlock::getSurfaceGeometricId(int i, int j, int k, int face)
{
  auto block = getGeometryBlock(i, j, k);
  return block->getGeometricEntity(LocalGeometricEntity(2, face));
}

int MeshGeneratorMultiBlock::getVolumeGeometricId(int i, int j, int k)
{
  auto block = getGeometryBlock(i, j, k);
  return block->getGeometricEntity(LocalGeometricEntity(3, 0));
}

// indices are in the "middle block" coordinate system: the middle block has
// indices (0, 0, 0)
std::shared_ptr<BlockGeometry> MeshGeneratorMultiBlock::getGeometryBlock(int i, int j, int k)
{
  assert(i >= -int(m_meshspec.numel_minusx.size()) && i <= int(m_meshspec.numel_plusx.size()));
  assert(j >= -int(m_meshspec.numel_minusy.size()) && j <= int(m_meshspec.numel_plusy.size()));
  assert(k >= -int(m_meshspec.numel_minusz.size()) && k <= int(m_meshspec.numel_plusz.size()));

  int iprime = i + m_meshspec.numel_minusx.size();
  int jprime = j + m_meshspec.numel_minusy.size();
  int kprime = k + m_meshspec.numel_minusz.size();

  return m_geometry_blocks[iprime][jprime][kprime];
}

void MeshGeneratorMultiBlock::reformatBlockData()
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

  m_middle_block_indices[0] = m_meshspec.numel_minusx.size();
  m_middle_block_indices[1] = m_meshspec.numel_minusy.size();
  m_middle_block_indices[2] = m_meshspec.numel_minusz.size();
}

void MeshGeneratorMultiBlock::computeBlockMeshSpecs()
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
        //std::cout << "block " << i << ", " << j << ", " << k << " meshspec = " << spec_ijk << std::endl;
        //std::cout << std::boolalpha << "create_block = " << m_meshspec.create_blocks[i][j][k] << std::endl;

        current_lower_corner_coords[2] += m_thicknesses[2][k];
      }

      current_lower_corner_coords[1] += m_thicknesses[1][j];
    }
    current_lower_corner_coords[0] += m_thicknesses[0][i];
  }
}

std::array<Real, 3> MeshGeneratorMultiBlock::getLowerCornerCoords()
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

void MeshGeneratorMultiBlock::initializeGeometry()
{
  gmi_register_null();
  m_gmi_topo = std::make_shared<mesh_gmi::GMITopo>();
  createGeometryBlocks();
  m_gmi_model = mesh_gmi::createGMITopo(m_gmi_topo);
  //m_gmi_model = gmi_load(".null");
}

void MeshGeneratorMultiBlock::createGeometryBlocks()
{
  m_geometry_blocks.resize(boost::extents[m_numels[0].size()][m_numels[1].size()][m_numels[2].size()]);
  for (int i=0; i < m_numels[0].size(); ++i)
    for (int j=0; j < m_numels[1].size(); ++j)
      for (int k=0; k < m_numels[2].size(); ++k)
      {
        if (!m_meshspec.create_blocks[i][j][k])
          continue;

        //std::cout << "\ncreating geometric block " << i << ", " << j << ", " << k << std::endl;
        auto surrounding_blocks = getSurroundingGeometryBlocks(i, j, k);
        m_geometry_blocks[i][j][k] = std::make_shared<BlockGeometry>(surrounding_blocks, m_geometric_id_gen, m_gmi_topo);
      }
}


void MeshGeneratorMultiBlock::createMeshBlocks()
{
  m_mesh_blocks.resize(boost::extents[m_numels[0].size()][m_numels[1].size()][m_numels[2].size()]);
  for (int i=0; i < m_numels[0].size(); ++i)
    for (int j=0; j < m_numels[1].size(); ++j)
      for (int k=0; k < m_numels[2].size(); ++k)
      {
        if (!m_meshspec.create_blocks[i][j][k])
          continue;

        //std::cout << "\ncreating mesh block " << i << ", " << j << ", " << k << std::endl;

        auto surrounding_blocks = getSurroundingMeshBlocks(i, j, k);
        m_mesh_blocks[i][j][k] = std::make_shared<MeshBlock>(m_mesh, m_meshspecs[i][j][k],
                                    m_geometry_blocks[i][j][k], surrounding_blocks);
        m_mesh_blocks[i][j][k]->generateMeshBlock();
      }
}

std::vector<std::shared_ptr<BlockGeometry>> MeshGeneratorMultiBlock::getSurroundingGeometryBlocks(int i, int j, int k)
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

std::vector<std::shared_ptr<MeshBlock>> MeshGeneratorMultiBlock::getSurroundingMeshBlocks(int i, int j, int k)
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

bool MeshGeneratorMultiBlock::isBlockValid(int i, int j, int k)
{
    return i >= 0 && i < m_numels[0].size() &&
            j >= 0 && j < m_numels[1].size() &&
            k >= 0 && k < m_numels[2].size() &&
            m_meshspec.create_blocks[i][j][k];
}


}  // namespace