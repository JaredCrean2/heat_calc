#include "mesh/mesh_block.h"

namespace Mesh {


void MeshBlock::generateMeshBlock()
{
  createVertices();
  createEdges();
  createFaces();
  createElements();
}


void MeshBlock::createVertices()
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


apf::MeshEntity* MeshBlock::getMeshVert(LocalGeometricEntity local_ge, int i, int j, int k)
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


LocalGeometricEntity MeshBlock::getVertModelEntity(const int i, const int j, const int k)
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


Point MeshBlock::computeVertCoords(const int i, const int j, const int k)
{
  double dx = (m_meshspec.xmax - m_meshspec.xmin)/m_meshspec.nx;
  double dy = (m_meshspec.ymax - m_meshspec.ymin)/m_meshspec.ny;
  double dz = (m_meshspec.zmax - m_meshspec.zmin)/m_meshspec.nz;

  Point p{i*dx + m_coords_min[0], j*dy + m_coords_min[1], k*dz + m_coords_min[2]};
  return p; //TODO: mapping function
  //return m_func(p);
}


void MeshBlock::createEdges()
{
  //std::cout << "\nMaking x direction edges" << std::endl;
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

  //std::cout << "\nMaking y direction edges" << std::endl;
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

  //std::cout << "\nMaking z direction edges" << std::endl;
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
LocalGeometricEntity MeshBlock::getEdgeModelEntity_x(const int i, const int j, const int k)
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
LocalGeometricEntity MeshBlock::getEdgeModelEntity_y(const int i, const int j, const int k)
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
LocalGeometricEntity MeshBlock::getEdgeModelEntity_z(const int i, const int j, const int k)
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

apf::MeshEntity* MeshBlock::getCommonEdge(apf::MeshEntity* vert1, apf::MeshEntity* vert2)
{
    for (int j=0; j < m_mesh->countUpward(vert1); ++j)
      for (int k=0; k < m_mesh->countUpward(vert2); ++k)
        if (m_mesh->getUpward(vert1, j) == m_mesh->getUpward(vert2, k))
          return m_mesh->getUpward(vert1, j);

    throw std::runtime_error("unable to find common edge");
}

void MeshBlock::printEdges(apf::MeshEntity* verts[4])
{
  for (int i=0; i < 4; ++i)
  {
    int i2 = (i + 1) % 4;
    apf::MeshEntity* edge = getCommonEdge(verts[i], verts[i2]);
    apf::ModelEntity* me = m_mesh->toModel(edge);
    std::cout << "edge " << i << " = " << edge << ", classified on geometric entity " << m_mesh->getModelType(me) << ", " << m_mesh->getModelTag(me) << std::endl;
  }
}

void MeshBlock::printVerts(apf::MeshEntity* verts[4])
{
  std::cout << "verts = " << std::endl;
  for (int i=0; i < 4; ++i)
    std::cout << verts[i] << ", ";
  std::cout << std::endl;
}

void MeshBlock::createFaces()
{
  apf::MeshEntity* verts[4];
  //std::cout << "\ncreating xz faces" << std::endl;
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

  //std::cout << "\ncreating yz faces" << std::endl;
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

  //std::cout << "\ncreating xy faces" << std::endl;
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
LocalGeometricEntity MeshBlock::getFaceModelEntity_xz(const int i, const int j, const int k)
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
LocalGeometricEntity MeshBlock::getFaceModelEntity_yz(const int i, const int j, const int k)
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
LocalGeometricEntity MeshBlock::getFaceModelEntity_xy(const int i, const int j, const int k)
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

void MeshBlock::createElements()
{
  //std::cout << "\nCreating elements" << std::endl;
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
}  // namespace