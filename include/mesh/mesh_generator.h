// generates structured meshes on a mapped cube

#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include <iostream>
#include "ProjectDefs.h"
#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "gmi_null.h"
#include "mesh/mesh_input.h"

namespace Mesh
{

Point identity(const Point& point);

// creates a mesh that can be mapped to a cube
// face geometry classification:
//    0 = yminus
//    1 = xplus
//    2 = yplus
//    3 = xminus
//    4 = zminus
//    5 = zplus
template <typename T>
class MeshGenerator
{
  public:
    // func must be Point = func(const Point&)
    MeshGenerator(const MeshSpec& mesh_spec, T func) :
      m_func(func),
      m_meshspec(mesh_spec),
      m_verts(boost::extents[mesh_spec.nx+1][mesh_spec.ny+1][mesh_spec.nz+1])
  {}

    apf::Mesh2* generate()
    {
      gmi_register_null();
      m_gmi_model = gmi_load(".null");
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
      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {
            Point p = computeVertCoords(i, j, k);
            apf::Vector3 p2(p.x, p.y, p.z);
            auto me = getVertModelEntity(i, j, k);
            m_verts[i][j][k] = m_mesh->createVert(me);
            m_mesh->setPoint(m_verts[i][j][k], 0, p2);
          }
    }

    apf::ModelEntity* getVertModelEntity(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;
      int ymax = m_meshspec.ny;
      int zmax = m_meshspec.nz;

      ModelEntitySpec me;
      // findModelEntity(dim, tag)
      // z = zmin vertices
      if (i == 0 && j == 0 && k == 0)
        me = ModelEntitySpec(0, 0);
      else if (i == xmax && j == 0 && k == 0)
        me = ModelEntitySpec(0, 1);
      else if (i == xmax && j == ymax && k == 0)
        me = ModelEntitySpec(0, 2);
      else if (i == 0 && j == ymax && k == 0)
        me = ModelEntitySpec(0, 3);
      // z = zmax vertices
      else if (i == 0 && j == 0 && k == zmax)
        me = ModelEntitySpec(0, 4);
      else if (i == xmax && j == 0 && k == zmax)
        me = ModelEntitySpec(0, 5);
      else if (i == xmax && j == ymax && k == zmax)
        me = ModelEntitySpec(0, 6);
      else if (i == 0 && j == ymax && k == zmax)
        me = ModelEntitySpec(0, 7);
      // edges on z = zmin quad
      else if (j == 0 && k == 0)
        me = ModelEntitySpec(1, 0);
      else if (i == xmax && k == 0)
        me = ModelEntitySpec(1, 1);
      else if (j == ymax && k == 0)
        me = ModelEntitySpec(1, 2);
      else if (i == 0 && k == 0)
        me = ModelEntitySpec(1, 3);
      // edges on z = zmax quad
      else if (j == 0 && k == zmax)
        me = ModelEntitySpec(1, 4);
      else if (i == xmax && k == zmax)
        me = ModelEntitySpec(1, 5);
      else if (j == ymax && k == zmax)
        me = ModelEntitySpec(1, 6);
      else if (i == 0 && k == zmax)
        me = ModelEntitySpec(1, 7);
      // edges parallel to z axis
      else if (i == 0 && j == 0)
        me = ModelEntitySpec(1, 8);
      else if (i == xmax && j == 0)
        me = ModelEntitySpec(1, 9);
      else if (i == xmax && j == ymax)
        me = ModelEntitySpec(1, 10);
      else if (i == 0 && j == ymax)
        me = ModelEntitySpec(1, 11);
      // faces
      else if (j == 0)
        me = ModelEntitySpec(2, 0);
      else if (i == xmax)
        me = ModelEntitySpec(2, 1);
      else if (j == ymax)
        me = ModelEntitySpec(2, 2);
      else if (i == 0)
        me = ModelEntitySpec(2, 3);
      else if (k == 0)
        me = ModelEntitySpec(2, 4);
      else if (k == zmax)
        me = ModelEntitySpec(2, 5);
      // interior
      else
        me = ModelEntitySpec(3, 0);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    Point computeVertCoords(const int i, const int j, const int k)
    {
      double dx = (m_meshspec.xmax - m_meshspec.xmin)/m_meshspec.nx;
      double dy = (m_meshspec.ymax - m_meshspec.ymin)/m_meshspec.ny;
      double dz = (m_meshspec.zmax - m_meshspec.zmin)/m_meshspec.nz;

      Point p{i*dx + m_meshspec.xmin, j*dy + m_meshspec.ymin, k*dz + m_meshspec.zmin};
      return m_func(p);
    }

    void createEdges()
    {
      apf::MeshEntity* verts[2];
      apf::ModelEntity* me;
      for (int i=0; i < m_meshspec.nx; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {
            // x direction
            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i+1][j][k];
            me = getEdgeModelEntity_x(i, j, k);
            m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
          }

      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {
            // y direction
            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i][j+1][k];
            me = getEdgeModelEntity_y(i, j, k);
            m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
          }

      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz; ++k)
          {
            // z direction
            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i][j][k+1];
            me = getEdgeModelEntity_z(i, j, k);
            m_mesh->createEntity(apf::Mesh::EDGE, me, verts);
          }
    }

    // gets model entity for edges parallel to the x axis
    apf::ModelEntity* getEdgeModelEntity_x(const int i, const int j, const int k)
    {
      //int xmax = m_meshspec.nx + 1;
      int ymax = m_meshspec.ny;
      int zmax = m_meshspec.nz;

      ModelEntitySpec me;
      // geometric edges
      if (j == 0 && k == 0)
        me = ModelEntitySpec(1, 0);
      else if (j == ymax && k == 0)
        me = ModelEntitySpec(1, 2);
      else if (j == 0 && k == zmax)
        me = ModelEntitySpec(1, 4);
      else if (j == ymax && k == zmax)
        me = ModelEntitySpec(1, 6);
      // faces
      else if (j == 0)
        me = ModelEntitySpec(2, 0);
      else if (k == 0)
        me = ModelEntitySpec(2, 4);
      else if (j == ymax)
        me = ModelEntitySpec(2, 2);
      else if (k == zmax)
        me = ModelEntitySpec(2, 5);
      // interior
      else
        me = ModelEntitySpec(3, 0);


      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for edges parallel to the y axis
    apf::ModelEntity* getEdgeModelEntity_y(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;
      //int ymax = m_meshspec.ny;
      int zmax = m_meshspec.nz;

      ModelEntitySpec me;
      // geometric edges
      if (i == 0 && k == 0)
        me = ModelEntitySpec(1, 3);
      else if (i == xmax && k == 0)
        me = ModelEntitySpec(1, 1);
      else if (i == 0 && k == zmax)
        me = ModelEntitySpec(1, 7);
      else if (i == xmax && k == zmax)
        me = ModelEntitySpec(1, 5);
      // faces
      else if (i == 0)
        me = ModelEntitySpec(2, 3);
      else if (k == 0)
        me = ModelEntitySpec(2, 4);
      else if (i == xmax)
        me = ModelEntitySpec(2, 1);
      else if (k == zmax)
        me = ModelEntitySpec(2, 5);
      // interior
      else
        me = ModelEntitySpec(3, 0);
 
      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for edges parallel to the z axis
    apf::ModelEntity* getEdgeModelEntity_z(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;
      int ymax = m_meshspec.ny;
      //int zmax = m_meshspec.nz;

      ModelEntitySpec me;
      // geometric edges
      if (j == 0 && i == 0)
        me = ModelEntitySpec(1, 8);
      else if (j == ymax && i == 0)
        me = ModelEntitySpec(1, 11);
      else if (j == 0 && i == xmax)
        me = ModelEntitySpec(1, 9);
      else if (j == ymax && i == xmax)
        me = ModelEntitySpec(1, 10);
      // faces
      else if (j == 0)
        me = ModelEntitySpec(2, 0);
      else if (i == 0)
        me = ModelEntitySpec(2, 3);
      else if (j == ymax)
        me = ModelEntitySpec(2, 2);
      else if (i == xmax)
        me = ModelEntitySpec(2, 1);
      // interior
      else
        me = ModelEntitySpec(3, 0);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    void createFaces()
    {
      apf::MeshEntity* verts[4];
      apf::ModelEntity* me;
      for (int i=0; i < m_meshspec.nx; ++i)
        for (int j=0; j < m_meshspec.ny+1; ++j)
          for (int k=0; k < m_meshspec.nz; ++k)
          {
            // xz plane
            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i+1][j][k];
            verts[2] = m_verts[i+1][j][k+1];
            verts[3] = m_verts[i][j][k+1];
            me = getFaceModelEntity_xz(i, j, k);
            apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
          }

      for (int i=0; i < m_meshspec.nx+1; ++i)
        for (int j=0; j < m_meshspec.ny; ++j)
          for (int k=0; k < m_meshspec.nz; ++k)
          {
            // yz plane
            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i][j][k+1];
            verts[2] = m_verts[i][j+1][k+1];
            verts[3] = m_verts[i][j+1][k];
            me = getFaceModelEntity_yz(i, j, k);
            apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
          }

      for (int i=0; i < m_meshspec.nx; ++i)
        for (int j=0; j < m_meshspec.ny; ++j)
          for (int k=0; k < m_meshspec.nz+1; ++k)
          {
            // xy plane
            verts[0] = m_verts[i][j][k];
            verts[1] = m_verts[i+1][j][k];
            verts[2] = m_verts[i+1][j+1][k];
            verts[3] = m_verts[i][j+1][k];
            me = getFaceModelEntity_xy(i, j, k);
            apf::buildElement(m_mesh, me, apf::Mesh::QUAD, verts);
          }

    }

    // gets model entity for faces parallel to the xz plane
    apf::ModelEntity* getFaceModelEntity_xz(const int i, const int j, const int k)
    {
      int ymax = m_meshspec.ny;

      ModelEntitySpec me;
      if (j == 0)
        me = ModelEntitySpec(2, 0);
      else if (j == ymax)
        me = ModelEntitySpec(2, 2);
      else
        me = ModelEntitySpec(3, 0);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for faces parallel to the yz plane
    apf::ModelEntity* getFaceModelEntity_yz(const int i, const int j, const int k)
    {
      int xmax = m_meshspec.nx;

      ModelEntitySpec me;
      if (i == 0)
        me = ModelEntitySpec(2, 3);
      else if (i == xmax)
        me = ModelEntitySpec(2, 1);
      else
        me = ModelEntitySpec(3, 0);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    // gets model entity for faces parallel to the xy plane
    apf::ModelEntity* getFaceModelEntity_xy(const int i, const int j, const int k)
    {
      int zmax = m_meshspec.nz;

      ModelEntitySpec me;
      if (k == 0)
        me = ModelEntitySpec(2, 4);
      else if (k == zmax)
        me = ModelEntitySpec(2, 5);
      else
        me = ModelEntitySpec(3, 0);

      return m_mesh->findModelEntity(me.dim, me.tag);
    }

    void createElements()
    {
      apf::MeshEntity* verts[8];
      apf::ModelEntity* me = m_mesh->findModelEntity(3, 0);
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

    T m_func; // mapping function
    MeshSpec m_meshspec;
    ArrayType<apf::MeshEntity*, 3> m_verts;
    gmi_model* m_gmi_model;
    apf::Mesh2* m_mesh;
};

template <typename T>
MeshGenerator<typename std::decay<T>::type>
make_mesh_generator(const MeshSpec& meshspec, T&& func=&identity)
{
  return {meshspec, std::forward<T>(func)};
}

} // namespace

#endif  // header guard
