#include "gtest/gtest.h"

#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC_defs.h"
#include "discretization/disc_vector.h"
#include "discretization/source_term.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/basis_vals.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern_mesh.h"
#include "mesh_helper.h"

#include "mpi.h"
#include "petscerror.h"


namespace {

  class HeatMMSConvergenceTester : public StandardDiscSetup,
                                   public testing::Test
  {
    protected:
      using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

      HeatMMSConvergenceTester()
      {
        setup();
      }

      using StandardDiscSetup::setup;

      virtual void setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& spec,
                         const linear_system::LargeMatrixOptsPetsc& matrix_opts,
                         const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true})
      {
        StandardDiscSetup::setup(quad_degree, sol_degree, spec, is_surf_dirichlet);

        heat        = std::make_shared<Heat::HeatEquation>(disc);
        u_vec       = makeDiscVector(disc);
        u_solve_vec = makeDiscVector(disc);
        res_vec     = makeDiscVector(disc);

        auto num_dofs     = disc->getDofNumbering()->getNumDofs();
        std::cout << "num_dofs = " << num_dofs << std::endl;
        this->matrix_opts = matrix_opts;
        auto sparsity     = std::make_shared<linear_system::SparsityPatternMesh>(mesh);
        mat               = std::make_shared<linear_system::LargeMatrixPetsc>(num_dofs, num_dofs, matrix_opts, sparsity);
        assembler         = std::make_shared<linear_system::Assembler>(disc, mat);
      }

      template <typename Tex, typename Tderiv, typename Tsrc>
      void setSolution(Tex ex_sol, Tderiv deriv, Tsrc src)
      {
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

      void computeErrorNorm(DiscVectorPtr u_ex, DiscVectorPtr u_sol)
      {
        int num_dofs = u_ex->getNumDofs();
        auto err = makeDiscVector(disc);
        auto tmp = makeDiscVector(disc);

        auto& u_ex_vec = u_ex->getVector();
        auto& u_sol_vec = u_sol->getVector();
        auto& err_vec = err->getVector();
        for (int i=0; i < num_dofs; ++i)
        {
          err_vec[i] = u_ex_vec[i] - u_sol_vec[i];          
        }
        err->markVectorModified();

        heat->applyMassMatrix(err, tmp);
        if (!tmp->isVectorCurrent())
          tmp->syncArrayToVector();
        
        Real val_l2 = 0, val_max = 0;
        auto& tmp_vec = tmp->getVector();
        for (int i=0; i < num_dofs; ++i)
        {
          val_l2 += err_vec[i] * tmp_vec[i];
          val_max = std::max(val_max, std::abs(err_vec[i]));
        }
        val_l2 = std::sqrt(val_l2);

        recordError(val_l2, val_max);
      }

      void recordError(Real err_l2, Real err_max, Real h)
      {
        std::cout << "recording err_l2 = " << err_l2 << ", err_max = " << err_max << " with h = " << h << std::endl;
        m_errors_l2.push_back(err_l2);
        m_errors_max.push_back(err_max);
        m_h_values.push_back(h);
      }

      void recordError(Real err_l2, Real err_max)
      {
        int nelems = spec.nx * spec.ny * spec.nz;
        Real h = std::pow(1.0/nelems, 1.0/3.0);
        recordError(err_l2, err_max, h);
      }

      ArrayType<Real, 2> computeSlopes()
      {
        assert(m_errors_l2.size() >= 2);
        ArrayType<Real, 2> slopes(boost::extents[m_errors_l2.size()][2]);
        slopes[0][0] = 0; slopes[0][1] = 0;
        for (unsigned int i=1; i < m_errors_l2.size(); ++i)
        {
          slopes[i][0] = std::log(m_errors_l2[i]/m_errors_l2[i-1]) / std::log(m_h_values[i]/m_h_values[i-1]);
          slopes[i][1] = std::log(m_errors_max[i]/m_errors_max[i-1]) / std::log(m_h_values[i]/m_h_values[i-1]);
        }

        return slopes;
      }

      std::vector<Real> getErrorsL2() const { return m_errors_l2; }
      std::vector<Real> getErrorsMax() const { return m_errors_max; }

      std::vector<Real> getHValues() const { return m_h_values; }

      HeatPtr heat;
      DiscVectorPtr u_vec;
      DiscVectorPtr u_solve_vec;
      DiscVectorPtr res_vec;
      linear_system::LargeMatrixOptsPetsc matrix_opts;
      linear_system::LargeMatrixPtr mat;
      linear_system::AssemblerPtr assembler;

    private:
      std::vector<Real> m_errors_l2;
      std::vector<Real> m_errors_max;
      std::vector<Real> m_h_values;
  };


  linear_system::LargeMatrixOptsPetsc get_options()
  {
    linear_system::LargeMatrixOptsPetsc opts;
    opts.is_structurally_symmetric = false;
    opts.is_value_symmetric        = false;
    opts.factor_in_place           = false;
    opts.petsc_opts["ksp_atol"] = "1e-15";
    opts.petsc_opts["ksp_rtol"] = "1e-50";
    opts.petsc_opts["ksp_monitor"] = "";
    //opts.petsc_opts["info"] = "";
    //PetscOptionsSetValue(nullptr, "-log_view", "");
    //PetscOptionsSetValue(nullptr, "-on_error_abort", "");
    //PetscPushErrorHandler(PetscAbortErrorHandler, nullptr);

    return opts;
  }
}

TEST(VerificaitonTest, Foo)
{
  EXPECT_TRUE(true);
}

TEST_F(HeatMMSConvergenceTester, Exponential)
{
  const int sol_degree = 1;
  auto opts = get_options();
  auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                      { return std::exp(x + y + z); };

  auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> std::array<Real, 3>
                      { return {std::exp(x + y + z), std::exp(x + y + z), std::exp(x + y + z)}; };
  auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { return -3*std::exp(x + y + z); };

  int nmeshes = 3;
  int nelem = 3;
  for (int i=0; i < nmeshes; ++i)
  {
    nelem = (i + 1) * 3;
    auto meshspec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, nelem, nelem, nelem);

    std::cout << "mesh " << i << std::endl;
    setup(2*sol_degree, sol_degree, meshspec, opts);
    setSolution(ex_sol_l, deriv_l, src_func_l);
    int num_dofs = u_vec->getNumDofs();

    heat->computeRhs(u_vec, 0.0, res_vec);
    res_vec->syncArrayToVector();

    heat->computeJacobian(u_vec, 0.0, assembler);

    mat->finishMatrixAssembly();
    mat->factor();
    solve(mat, res_vec, u_solve_vec);

    auto& u_solve = u_solve_vec->getVector();
    for (int j=0; j < num_dofs; ++j)
      u_solve[j] = u_vec->getVector()[j] - u_solve[j];
    u_solve_vec->markVectorModified();
    

    computeErrorNorm(u_vec, u_solve_vec);
  }

  auto slopes = computeSlopes();
  auto errors_l2 = getErrorsL2();
  auto errors_max = getErrorsMax();

  for (unsigned int i=0; i < slopes.size(); ++i)
  {
    Real ratio_l2 = 0, ratio_max = 0;
    if (i > 0)
    {
      ratio_l2 = errors_l2[i-1] / errors_l2[i];
      ratio_max = errors_max[i-1] / errors_max[i];
    }

    std::cout << "mesh " << i << ": slope_l2 = " << slopes[i][0] << ", slope_max = " << slopes[i][1] << ", ratio_l2 = " << ratio_l2 << ", ratio_max = " << ratio_max << std::endl;
  }

  //TODO: make mesh finer
}