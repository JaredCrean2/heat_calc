#include "gtest/gtest.h"

#include "mesh_helper.h"
#include "discretization/surface_discretization.h"
#include "mesh/mesh_generator.h"

TEST(SurfDisc, Normals)
{
  auto spec = getStandardMeshSpec();
  spec.nx = 4; spec.ny = 5, spec.nz = 6;
  spec.xmax = 1; spec.ymax = 1; spec.zmax = 1;
  auto mesh = makeStandardMesh(spec);

  Quadrature quad = getGaussianQuadrature(1);
  Mesh::VolumeGroup& vol_group = mesh->getElements(0);
  auto vol_disc = std::make_shared<VolumeDiscretization>(vol_group, quad);
  std::vector<std::shared_ptr<VolumeDiscretization>> vol_discs{vol_disc};

  std::vector<std::shared_ptr<SurfaceDiscretization>> surf_discs;
  for (Index i=0; i < mesh->getNumSurfaces(); ++i)
  {
    auto surf_i = std::make_shared<SurfaceDiscretization>(mesh->getFaces(i), quad, vol_discs); 
    surf_discs.push_back(surf_i);
  }


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


