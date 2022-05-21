#include "gtest/gtest.h"
#include "test_helper.h"
//#include "mesh/mesh.h"
#include "mesh_helper.h"
#include "discretization/surface_discretization.h"
#include "mesh/mesh_generator.h"
#include "quadrature.h"

namespace {
class SurfaceTester : public ::testing::Test,
                      public StandardMeshSetup

{
  protected:
    SurfaceTester()
    {
      SERIAL_ONLY_RETURN();
      setup();
    }
};

}

TEST_F(SurfaceTester, Counts)
{
  SERIAL_ONLY();

  int sol_degree = 1;
  int quad_degree = 3;
  setup(quad_degree, sol_degree);
  for (auto& surf : surf_discs)
  {
    EXPECT_EQ(surf->getNumCoordPtsPerFace(), 4);
    EXPECT_EQ(surf->getNumQuadPtsPerFace(), 4);
    EXPECT_EQ(surf->getNumSolPtsPerFace(), 4);
  }

  quad_degree = 5;
  setup(quad_degree, sol_degree);
  for (auto& surf : surf_discs)
  {
    EXPECT_EQ(surf->getNumCoordPtsPerFace(), 4);
    EXPECT_EQ(surf->getNumQuadPtsPerFace(), 9);
    EXPECT_EQ(surf->getNumSolPtsPerFace(), 4);
  }


  sol_degree = 2;
  quad_degree = 5;
  setup(quad_degree, sol_degree);
  for (auto& surf : surf_discs)
  {
    EXPECT_EQ(surf->getNumCoordPtsPerFace(), 4);
    EXPECT_EQ(surf->getNumQuadPtsPerFace(), 9);
    EXPECT_EQ(surf->getNumSolPtsPerFace(), 9);
  }


  sol_degree = 3;
  quad_degree = 5;
  setup(quad_degree, sol_degree);
  for (auto& surf : surf_discs)
  {
    EXPECT_EQ(surf->getNumCoordPtsPerFace(), 4);
    EXPECT_EQ(surf->getNumQuadPtsPerFace(), 9);
    EXPECT_EQ(surf->getNumSolPtsPerFace(), 16);
  }
}


TEST_F(SurfaceTester, Normals)
{
  SERIAL_ONLY();

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

  for (unsigned int i=0; i < surf_discs.size(); ++i)
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
  SERIAL_ONLY();

  // the coordinates that vary over each geometric face
  ArrayType<int, 2> active_coords(boost::extents[6][2]);
  active_coords[0][0] = 0; active_coords[0][1] = 2;
  active_coords[1][0] = 1; active_coords[1][1] = 2;
  active_coords[2][0] = 0; active_coords[2][1] = 2;
  active_coords[3][0] = 1; active_coords[3][1] = 2;
  active_coords[4][0] = 0; active_coords[4][1] = 1;
  active_coords[5][0] = 0; active_coords[5][1] = 1;

  for (unsigned int i=0; i < surf_discs.size(); ++i)
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

        // compute exact result
        int other_active_dim = (active_coords[i][d] + 1) % 2;
        double area = mesh_dim_maxs[other_active_dim] - mesh_dim_mins[other_active_dim];

        double val_exact = area*(poly_antideriv(mesh_dim_maxs[d], degree) - poly_antideriv(mesh_dim_mins[d], degree));
        EXPECT_FLOAT_EQ(val_sum, val_exact);
      }
  }
}

TEST_F(SurfaceTester, integrateFaceVector)
{
  SERIAL_ONLY();
  
  std::vector<int> active_coords{1, 0, 1, 0, 2, 2};
  std::vector<int> face_sign{-1, 1, 1, -1, -1, 1};

  for (unsigned int i=0; i < surf_discs.size(); ++i)
  {
    auto surf_i = surf_discs[i];
    auto& quad = surf_i->quad;
    auto active_coord = active_coords[i];  // normal vector component that is non-zero
    auto variable_coord = (active_coord + 1) % 3;  // coordinate that varies over the face
    auto other_variable_coord = (active_coord + 2) % 3;
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
        
      // compute exact result
      double area = mesh_dim_maxs[other_variable_coord] - mesh_dim_mins[other_variable_coord];
      double val_exact = area*(poly_antideriv(mesh_dim_maxs[variable_coord], degree) -
                                poly_antideriv(mesh_dim_mins[variable_coord], degree)) * face_sign[i];
      EXPECT_FLOAT_EQ(val_sum, val_exact);
    }
  }
}


