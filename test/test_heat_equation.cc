
#include "gtest/gtest.h"
#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC_defs.h"
#include "discretization/disc_vector.h"
#include "discretization/source_term.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/basis_vals.h"
#include "mesh_helper.h"

namespace {
  class HeatMMSTester : public StandardDiscSetup,
                        public testing::Test
  {
    protected:
      using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

      HeatMMSTester()
      {
        setup();
      }

      template <typename Tex, typename Tderiv, typename Tsrc>
      void setSolution(Tex ex_sol, Tderiv deriv, Tsrc src)
      {
        heat = std::make_shared<Heat::HeatEquation>(disc);
        u_vec = makeDiscVector(disc);
        res_vec = makeDiscVector(disc);


        for (int i=0; i < disc->getNumVolDiscs(); ++i)
        {
          heat->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
          heat->addVolumeGroupParams(Heat::VolumeGroupParams(1, 1, 1));
        }

        for (int i=0; i < disc->getNumSurfDiscs(); ++i)
        {
          auto surf = disc->getSurfDisc(i);
          if (surf->getIsDirichlet())
            heat->addDirichletBC(makeDirichletBCMMS(surf, ex_sol));
          else
            heat->addNeumannBC(makeNeumannBCMMS(surf, deriv));
        }


        res_vec->set(0);
        auto f = [&](Real x, Real y, Real z)
                    { return ex_sol(x, y, z, 0); };
        u_vec->setFunc(f);
      }

      HeatPtr heat;
      DiscVectorPtr u_vec;
      DiscVectorPtr res_vec;
    };

Real ex_sol(Real x, Real y, Real z, Real t, int degree)
{
  return std::pow(x, degree); // + std::pow(y, degree) + std::pow(z, degree);
}


Real ex_sol_deriv(Real x, Real y, Real z, Real t, int degree)
{
  if (degree == 0)
    return 0;
  else
    return degree * std::pow(x, degree - 1); //+ degree * std::pow(y, degree - 1) + degree * std::pow(z, degree-1);
}


Real src_func_dir(Real x, int degree)
{
  if (degree == 0 || degree == 1)
    return 0;
  else
    return -(degree * (degree - 1) * std::pow(x, degree-2));
}

Real src_func(Real x, Real y, Real z, Real t, int degree)
{
  return src_func_dir(x, degree); // + src_func_dir(y, degree) + src_func_dir(z, degree);
}

class StandardDiscTester : public StandardDiscSetup,
                           public ::testing::Test
{
};

void getQuadTPIndices(SurfDiscPtr surf, int quad_point, int& i, int& j)
{
  int num_quad_pts_per_dir = surf->quad_tp_nodemap.shape()[0];
  for (i=0; i < num_quad_pts_per_dir; ++i)
    for (j=0; j < num_quad_pts_per_dir; ++j)
      if (surf->quad_tp_nodemap[i][j] == quad_point)
        return;

  throw std::runtime_error("unable to find quadrature point");
}

}  // namespace


TEST_F(HeatMMSTester, Constant)
{
  auto ex_sol = [] (Real x, Real y, Real z, Real t)
                   { return 1; };

  auto deriv = [] (Real x, Real y, Real z, Real t)
                   { return 1; };

  auto src_func = [] (Real x, Real y, Real z, Real t)
                     { return 0; };
  
  //setup();
  setSolution(ex_sol, deriv, src_func);

  heat->computeRhs(u_vec, 0.0, res_vec);
  res_vec->syncArrayToVector();

  auto& vec = res_vec->getVector();
  for (int i=0; i < vec.shape()[0]; ++i)
  {
    EXPECT_LE(std::abs(vec[i]), 1e-13);
  }
}


TEST_F(HeatMMSTester, PolynomialExactnessDirichlet)
{
  for (int sol_degree=1; sol_degree <= 3; ++sol_degree)
  {
    std::cout << "testing sol degree " << sol_degree << std::endl;
    for (int degree=0; degree <= sol_degree; ++degree)
    {
      std::cout << "  testing polynomial degree " << degree << std::endl;
      auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return ex_sol(x, y, z, t, degree); };

      auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return ex_sol_deriv(x, y, z, t, degree); };
      auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                            { return src_func(x, y, z, t, degree); };

      setup(2*sol_degree, sol_degree);
      setSolution(ex_sol_l, deriv_l, src_func_l);

      heat->computeRhs(u_vec, 0.0, res_vec);
      res_vec->syncArrayToVector();

      auto& vec = res_vec->getVector();
      for (int i=0; i < vec.shape()[0]; ++i)
      {
        EXPECT_LE(std::abs(vec[i]), 1e-12);
      }
    }
  }
}


TEST_F(HeatMMSTester, PolynomialExactnessNeumann)
{
  std::vector<bool> dirichlet_surfs = {false, true, false, true, true, true};
  for (int sol_degree=1; sol_degree <= 1; ++sol_degree)
  {
    std::cout << "testing sol degree " << sol_degree << std::endl;
    //for (int degree=0; degree <= sol_degree; ++degree)
    for (int degree=0; degree <= 1; ++degree)
    {
      std::cout << "  testing polynomial degree " << degree << std::endl;
      auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return ex_sol(x, y, z, t, degree); };

      auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return ex_sol_deriv(x, y, z, t, degree); };

      auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                            { return src_func(x, y, z, t, degree); };

      setup(2*sol_degree, sol_degree, dirichlet_surfs);
      setSolution(ex_sol_l, deriv_l, src_func_l);

      heat->computeRhs(u_vec, 0.0, res_vec);
      res_vec->syncArrayToVector();

      auto& vec = res_vec->getVector();
      for (int i=0; i < vec.shape()[0]; ++i)
      {
        EXPECT_LE(std::abs(vec[i]), 1e-12);
      }
    }
  }
}


TEST_F(StandardDiscTester, BasisVals2D_interpolation)
{
  for (int sol_degree=1; sol_degree <= 3; ++sol_degree)
  {
    int quad_degree = 2*sol_degree;
    std::vector<bool> is_surf_dirichlet = {false, false, false, false, false, false};
    setup(quad_degree, sol_degree, is_surf_dirichlet);

    for (int degree=0; degree <= sol_degree; ++degree)
      for (int surf_i=0; surf_i < 6; ++surf_i)
      {
        //std::cout << "\ntesting surface " << surf_i << ", degree " << degree << std::endl;
        auto surf = disc->getSurfDisc(surf_i);
        auto& quad = surf->quad;
        Heat::BasisVals2D basis(surf->face_group.getTPMapperSol(), quad.getPoints(), surf->face_group.getFaceNodesSol(), surf->face_group.ref_el_sol);
        ArrayType<Real, 2> quad_coords(boost::extents[surf->getNumQuadPtsPerFace()][3]);
        ArrayType<Real, 2> sol_coords(boost::extents[surf->getNumSolPtsPerFace()][3]);
        ArrayType<Real, 1> sol_vals(boost::extents[surf->getNumSolPtsPerFace()]);

        for (int face=0; face < surf->getNumFaces(); ++face)
        {
          auto& face_spec = surf->getFaceSpec(face);
          surf->getFaceQuadCoords(face, quad_coords);
          surf->getFaceSolCoords(face, sol_coords);

          for (int k=0; k < surf->getNumSolPtsPerFace(); ++k)
            sol_vals[k] = ex_sol(sol_coords[k][0], sol_coords[k][1], sol_coords[k][2], 0.0, degree);

          // interpolate to quadrature points
          for (int q=0; q < surf->getNumQuadPtsPerFace(); ++q)
          {
            Real val = 0.0;
            for (int i=0; i < surf->getNumSolPtsPerFace(); ++i)
            {
              // get tp indices
              int qi, qj;
              getQuadTPIndices(surf, q, qi, qj);
              val += basis.getValue(face_spec.face, i, qi, qj) * sol_vals[i];
            }

            Real val_ex = ex_sol(quad_coords[q][0], quad_coords[q][1], quad_coords[q][2], 0.0, degree);
            //std::cout << "val = " << val << ", val_ex = " << val_ex << ", diff = " << val - val_ex << std::endl;

            EXPECT_NEAR(val, val_ex, 1e-13);

          }

        }
          
      }
  }
}
