#include "gtest/gtest.h"

#include "ProjectDefs.h"
#include "mesh_helper.h"
#include "discretization/volume_discretization.h"
#include "mesh/mesh_generator.h"

TEST(VolDisc, dxidx)
{
  auto spec = getStandardMeshSpec();
  spec.nx = 4; spec.ny = 5, spec.nz = 6;
  auto mesh = makeStandardMesh(spec);

  Quadrature quad = getGaussianQuadrature(1);
  Mesh::VolumeGroup& vol_group = mesh->getElements(0);
  VolumeDiscretization vol_disc(vol_group, quad);

  // xi is 0 to 1, the mesh is linear and axis aligned, so dxidx = 1/(xmax - xmin)/nx)
  
  double dxidx_ex[3] = {1.0/( (spec.xmax - spec.xmin)/spec.nx ),
                        1.0/( (spec.ymax - spec.ymin)/spec.ny ),
                        1.0/( (spec.zmax - spec.zmin)/spec.nz )};

  ArrayType<Real, 2> dxidx_ex_mat(boost::extents[3][3]);
  dxidx_ex_mat[0][0] = 0;            dxidx_ex_mat[0][1] = -dxidx_ex[1]; dxidx_ex_mat[0][2] = 0;
  dxidx_ex_mat[1][0] = dxidx_ex[0];  dxidx_ex_mat[1][1] = 0;            dxidx_ex_mat[1][2] = 0;
  dxidx_ex_mat[2][0] = 0;            dxidx_ex_mat[2][1] = 0;            dxidx_ex_mat[2][2] = dxidx_ex[2];

  

  for (int i=0; i < vol_disc.getNumElems(); ++i)
    for (int j=0; j < vol_disc.getNumQuadPtsPerElement(); ++j)
      for (int d1=0; d1 < 3; ++d1)
      {
        for (int d2=0; d2 < 3; ++d2)
        
        {
          EXPECT_NEAR(vol_disc.dxidx[i][j][d1][d2], dxidx_ex_mat[d1][d2], 1e-13);
          /*
          if (d1 == d2)
          {
            EXPECT_NEAR(vol_disc.dxidx[i][j][d1][d2], dxidx_ex[d1], 1e-13);
          } else
          {
            EXPECT_NEAR(vol_disc.dxidx[i][j][d1][d2], 0, 1e-13);
          }
          */
          
        }
      }
}

TEST(VolDisc, MatFuncs)
{
  ArrayType<Real, 2> A(boost::extents[3][3]), C(boost::extents[3][3]);

  A[0][0] = 1; A[1][0] = 4;  A[2][0] = 8;
  A[0][1] = 2; A[1][1] = 5;  A[2][1] = 8;
  A[0][2] = 3; A[1][2] = 6;  A[2][2] = 9;

  // determinant
  Real detA = computeDet3x3(A);
  EXPECT_FLOAT_EQ(detA, -3.0);

  // cofactor
  computeCofactor(A, C);
  ArrayType<Real, 2> C2(boost::extents[3][3]);
  C2[0][0] = -3; C2[0][1] =  12; C2[0][2] = -8;
  C2[1][0] =  6; C2[1][1] = -15; C2[1][2] =  8;
  C2[2][0] = -3; C2[2][1] =   6; C2[2][2] = -3;

  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      EXPECT_FLOAT_EQ(C[i][j], C2[i][j]);

  // inverse
  ArrayType<Real, 2> Ainv(boost::extents[3][3]);
  Ainv[0][0] =  1; Ainv[0][1] =  -2; Ainv[0][2] =  1;
  Ainv[1][0] = -4; Ainv[1][1] =   5; Ainv[1][2] = -2;
  Ainv[2][0] = 8.0/3.0; Ainv[2][1] =   -8.0/3.0; Ainv[2][2] = 1;

  computeInverse3x3(A, C);
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      EXPECT_FLOAT_EQ(Ainv[i][j], C[i][j]);
}


namespace {

class VolumeTester : public ::testing::Test,
                     public StandardMeshSetup

{
  protected:
    VolumeTester()
    {
      setup();
    }
};

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


TEST_F(VolumeTester, integrateVolumeScalar)
{
  ArrayType<Real, 2> quad_coords(boost::extents[vol_disc->getNumQuadPtsPerElement()][3]);
  ArrayType<Real, 1> vals(boost::extents[vol_disc->getNumQuadPtsPerElement()]);
  for (int degree=0; degree < 2*quad.getNumPoints() - 1; ++degree)
    for (int d=0; d < 3; ++d)
    {
      Real val = 0.0;
      for (int i=0; i < vol_disc->getNumElems(); ++i)
      {
        vol_disc->getVolumeQuadCoords(i, quad_coords);
        for (int j=0; j < vol_disc->getNumQuadPtsPerElement(); ++j)
          vals[j] = poly(quad_coords[j][d], degree);

        val += integrateVolumeScalar(vol_disc, i, vals);
      }

      // compute exact value
      Real area = 1;
      for (int d2=0; d2 < 3; ++d2)
        if (d2 != d)
          area *= mesh_dim_maxs[d2] - mesh_dim_mins[d2];

      Real val_exact = area*(poly_antideriv(mesh_dim_maxs[d], degree) - poly_antideriv(mesh_dim_mins[d], degree));
      EXPECT_FLOAT_EQ(val, val_exact);
    }
}

