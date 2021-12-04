#include "gtest/gtest.h"

#include "mesh/mesh.h"
#include "mesh_helper.h"
#include "discretization/surface_discretization.h"
#include "mesh/mesh_generator.h"
#include "quadrature.h"

namespace {
class SurfaceTester : public::testing::Test
{
  protected:
    SurfaceTester() :
      quad(getGaussianQuadrature(3))
    {
      spec = getStandardMeshSpec();
      //spec.nx = 4; spec.ny = 5, spec.nz = 6;
      spec.nx = 1; spec.ny = 1, spec.nz = 1;
      spec.xmax = 1; spec.ymax = 1; spec.zmax = 1;
      mesh = makeStandardMesh(spec);
 
      Mesh::VolumeGroup& vol_group = mesh->getElements(0);
      auto vol_disc = std::make_shared<VolumeDiscretization>(vol_group, quad);
      std::vector<std::shared_ptr<VolumeDiscretization>> vol_discs{vol_disc};

      //std::vector<std::shared_ptr<SurfaceDiscretization>> surf_discs;
      for (Index i=0; i < mesh->getNumSurfaces(); ++i)
      {
        auto surf_i = std::make_shared<SurfaceDiscretization>(mesh->getFaces(i), quad, vol_discs); 
        surf_discs.push_back(surf_i);
      }
    }

    Mesh::MeshSpec spec;
    std::shared_ptr<Mesh::MeshCG> mesh;
    Quadrature quad;
    std::vector<SurfDiscPtr> surf_discs;
};


}


TEST_F(SurfaceTester, Normals)
{
  // exact normals for the standard mesh
  ArrayType<Real, 2> normals_ex(boost::extents[6][3]);
  normals_ex[0][0] =  0; normals_ex[0][1] = -1; normals_ex[0][2] =  0;
  normals_ex[1][0] =  1; normals_ex[1][1] =  0; normals_ex[1][2] =  0;
  normals_ex[2][0] =  0; normals_ex[2][1] =  1; normals_ex[2][2] =  0;
  normals_ex[3][0] = -1; normals_ex[3][1] =  0; normals_ex[3][2] =  0;
  normals_ex[4][0] =  0; normals_ex[4][1] =  0; normals_ex[4][2] = -1;
  normals_ex[5][0] =  0; normals_ex[5][1] =  0; normals_ex[5][2] =  1;

  // exact face areas
  ArrayType<Real, 1> face_areas(boost::extents[6]);
  Real dx = (spec.xmax - spec.xmin)/spec.nx;
  Real dy = (spec.ymax - spec.ymin)/spec.ny;
  Real dz = (spec.zmax - spec.zmin)/spec.nz;
  face_areas[0] = dx*dz; face_areas[2] = dx*dz;
  face_areas[1] = dy*dz; face_areas[3] = dy*dz;
  face_areas[4] = dx*dy; face_areas[5] = dx*dy;

  for (Index i=0; i < mesh->getNumSurfaces(); ++i)
  {
    auto surf_i = surf_discs[i];

    ArrayType<Real, 1> ones(boost::extents[surf_i->getNumQuadPtsPerFace()]);
    for (int k=0; k < surf_i->getNumQuadPtsPerFace(); ++k)
      ones[k] = 1;

    for (int j=0; j < surf_i->getNumFaces(); ++j)
      for (int k=0; k < surf_i->getNumQuadPtsPerFace(); ++k)
      {
        // check face direction
        Real val = 0, mag1 = 0, mag_ex = 1;
        for (int d=0; d < 3; ++d)
        {
          val  += surf_i->normals[j][k][d] * normals_ex[i][d];
          mag1 += surf_i->normals[j][k][d] * surf_i->normals[j][k][d];
        }

        mag1 = std::sqrt(mag1);

        EXPECT_FLOAT_EQ(val, mag1 * mag_ex);

        // check magnitude (integrating normals gives face area)
        EXPECT_FLOAT_EQ(integrateFaceScalar(surf_i, j, ones), face_areas[i]);


      }
  }
}

namespace {
  double poly(const double x, const int n)
  {
    assert(n >= 0);
    return std::pow(x, n);
  }

  double poly_antideriv(const double x, const int n)
  {
    assert( n >= 0);      
    return std::pow(x, n+1)/(n+1);
  }
}

TEST_F(SurfaceTester, integrateFaceScalar)
{
  // the coordinates that vary over each geometric face
  ArrayType<int, 2> active_coords(boost::extents[6][2]);
  active_coords[0][0] = 0; active_coords[0][1] = 2;
  active_coords[1][0] = 1; active_coords[1][1] = 2;
  active_coords[2][0] = 0; active_coords[2][1] = 2;
  active_coords[3][0] = 1; active_coords[3][1] = 2;
  active_coords[4][0] = 0; active_coords[4][1] = 1;
  active_coords[5][0] = 0; active_coords[5][1] = 1;

  for (Index i=0; i < mesh->getNumSurfaces(); ++i)
  {
    auto surf_i = surf_discs[i];
    auto& quad = surf_i->quad;
    ArrayType<double, 2> face_coords(boost::extents[surf_i->getNumQuadPtsPerFace()][3]);
    ArrayType<double, 1> vals(boost::extents[surf_i->getNumQuadPtsPerFace()]);

    for (int degree=0; degree < 2*quad.getNumPoints() - 1; ++degree)
      for (int d=0; d < 2; ++d)
      {
        double val_sum = 0;
        for (int j=0; j < surf_i->getNumFaces(); ++j)
        {
          surf_i->getFaceQuadCoords(j, face_coords);
                      
          for (int k=0; k < surf_i->getNumQuadPtsPerFace(); ++k)
            vals[k] = poly(face_coords[k][active_coords[i][d]], degree);

          val_sum += integrateFaceScalar(surf_i, j, vals);
        }

        double val_exact = 1*(poly_antideriv(1.0, degree) - poly_antideriv(0.0, degree));
        EXPECT_FLOAT_EQ(val_sum, val_exact);
      }
  }
}

TEST_F(SurfaceTester, integrateFaceVector)
{
  std::vector<int> active_coords{1, 0, 1, 0, 2, 2};
  std::vector<int> face_sign{-1, 1, 1, -1, -1, 1};

  for (Index i=0; i < mesh->getNumSurfaces(); ++i)
  {
    auto surf_i = surf_discs[i];
    auto& quad = surf_i->quad;
    auto active_coord = active_coords[i];  // normal vector component that is non-zero
    auto variable_coord = (active_coord + 1) % 3;  // coordinate that varies over the face
    ArrayType<double, 2> face_coords(boost::extents[surf_i->getNumQuadPtsPerFace()][3]);
    ArrayType<double, 2> vals(boost::extents[surf_i->getNumQuadPtsPerFace()][3]);

    for (int k=0; k < surf_i->getNumQuadPtsPerFace(); ++k)
      for (int d=0; d < 3; ++d)
        if (d != active_coord)
          vals[k][d] = 666; // random data;

    for (int degree=0; degree < 2*quad.getNumPoints() - 1; ++degree)
    {
      double val_sum = 0.0;
      for (int j=0; j < surf_i->getNumFaces(); ++j)
      {
        surf_i->getFaceQuadCoords(j, face_coords);

        for (int k=0; k < surf_i->getNumQuadPtsPerFace(); ++k)
          vals[k][active_coord] = poly(face_coords[k][variable_coord], degree);

        val_sum += integrateFaceVector(surf_i, j, vals); 
      }

      double val_exact = 1*(poly_antideriv(1.0, degree) - poly_antideriv(0.0, degree)) * face_sign[i];
      EXPECT_FLOAT_EQ(val_sum, val_exact);
    }
  }
}


