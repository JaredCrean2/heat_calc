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
#include "time_solver/crank_nicolson.h"
#include "mesh_helper.h"

#include "mpi.h"
#include "petscerror.h"


namespace {

  class SlopeCalculator
  {
    public:
      void recordError(Real err_l2, Real err_max, Real h)
      {
        std::cout << "recording err_l2 = " << err_l2 << ", err_max = " << err_max << " with h = " << h << std::endl;
        m_errors_l2.push_back(err_l2);
        m_errors_max.push_back(err_max);
        m_h_values.push_back(h);
      }

      void recordError(Real err_l2, Real err_max, int nelems)
      {
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

      void resetErrors()
      {
        m_errors_l2.resize(0);
        m_errors_max.resize(0);
        m_h_values.resize(0);
      }

    private:
      std::vector<Real> m_errors_l2;
      std::vector<Real> m_errors_max;
      std::vector<Real> m_h_values;
  };

  class HeatMMSConvergenceTester : public StandardDiscSetup,
                                   public SlopeCalculator,
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

        auto num_dofs     = disc->getDofNumbering()->getNumOwnedDofs();
        this->matrix_opts = matrix_opts;
        auto sparsity     = std::make_shared<linear_system::SparsityPatternMesh>(mesh);
        mat               = std::make_shared<linear_system::LargeMatrixPetsc>(num_dofs, num_dofs, matrix_opts, sparsity);
        assembler         = std::make_shared<linear_system::Assembler>(disc, mat);
      }

      template <typename Tex, typename Tderiv, typename Tsrc>
      void setSolution(Tex ex_sol, Tderiv deriv, Tsrc src, const Heat::VolumeGroupParams& params = Heat::VolumeGroupParams{1, 1, 1})
      {
        for (int i=0; i < disc->getNumVolDiscs(); ++i)
        {
          heat->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
          heat->addVolumeGroupParams(params);
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
        std::vector<DofInt> owned_to_local_dofs;
        mesh->getOwnedLocalDofInfo(owned_to_local_dofs);
        int num_dofs = u_ex->getNumDofs();
        auto err     = makeDiscVector(disc);
        auto tmp     = makeDiscVector(disc);

        auto& u_ex_vec  = u_ex->getVector();
        auto& u_sol_vec = u_sol->getVector();
        auto& err_vec   = err->getVector();
        for (int i=0; i < num_dofs; ++i)
          err_vec[i] = u_ex_vec[i] - u_sol_vec[i];          

        err->markVectorModified();

        heat->applyMassMatrix(err, tmp);
        if (!tmp->isVectorCurrent())
          tmp->syncArrayToVector();
        
        double val_l2 = 0, val_max = 0;
        auto& tmp_vec = tmp->getVector();
        for (auto i : owned_to_local_dofs)
        {
          val_l2 += err_vec[i] * tmp_vec[i];
          val_max = std::max(val_max, std::abs(err_vec[i]));
        }

        double val_l2_global, val_max_global;
        MPI_Allreduce(&val_l2, &val_l2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&val_max, &val_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        val_l2_global = std::sqrt(val_l2_global);

        recordError(val_l2_global, val_max_global, spec.nx * spec.ny * spec.nz);
      }

      HeatPtr heat;
      DiscVectorPtr u_vec;
      DiscVectorPtr u_solve_vec;
      DiscVectorPtr res_vec;
      linear_system::LargeMatrixOptsPetsc matrix_opts;
      linear_system::LargeMatrixPtr mat;
      linear_system::AssemblerPtr assembler;
  };


  class UnsteadyHeatMMSConvergenceTester : public StandardDiscSetup,
                                           public SlopeCalculator,
                                           public testing::Test
  {
    protected:
      using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

      UnsteadyHeatMMSConvergenceTester()
      {
        setup();
      }

      using StandardDiscSetup::setup;

      virtual void setup(const int quad_degree, int sol_degree, const Mesh::MeshSpec& spec,
                         const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true})
      {
        StandardDiscSetup::setup(quad_degree, sol_degree, spec, is_surf_dirichlet);

        heat        = std::make_shared<Heat::HeatEquation>(disc);
        u_vec       = makeDiscVector(disc);
      }

      template <typename Tex, typename Tderiv, typename Tsrc>
      void setSolution(Tex ex_sol, Tderiv deriv, Tsrc src, const Heat::VolumeGroupParams& params = Heat::VolumeGroupParams{1, 1, 1})
      {
        for (int i=0; i < disc->getNumVolDiscs(); ++i)
        {
          heat->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
          heat->addVolumeGroupParams(params);
        }

        for (int i=0; i < disc->getNumSurfDiscs(); ++i)
        {
          auto surf = disc->getSurfDisc(i);
          if (surf->getIsDirichlet())
            heat->addDirichletBC(makeDirichletBCMMS(surf, ex_sol));
          else
            heat->addNeumannBC(makeNeumannBCMMS(surf, deriv));
        }

        auto f = [&](Real x, Real y, Real z)
                    { return ex_sol(x, y, z, 0); };
        u_vec->setFunc(f);
      }

      void computeErrorNorm(DiscVectorPtr u_ex, DiscVectorPtr u_sol)
      {
        std::vector<DofInt> owned_to_local_dofs;
        mesh->getOwnedLocalDofInfo(owned_to_local_dofs);
        int num_dofs = u_ex->getNumDofs();
        auto err     = makeDiscVector(disc);
        auto tmp     = makeDiscVector(disc);

        auto& u_ex_vec  = u_ex->getVector();
        auto& u_sol_vec = u_sol->getVector();
        auto& err_vec   = err->getVector();
        for (int i=0; i < num_dofs; ++i)
          err_vec[i] = u_ex_vec[i] - u_sol_vec[i];          

        err->markVectorModified();

        heat->applyMassMatrix(err, tmp);
        if (!tmp->isVectorCurrent())
          tmp->syncArrayToVector();
        
        double val_l2 = 0, val_max = 0;
        auto& tmp_vec = tmp->getVector();
        for (auto i : owned_to_local_dofs)
        {
          val_l2 += err_vec[i] * tmp_vec[i];
          val_max = std::max(val_max, std::abs(err_vec[i]));
        }

        double val_l2_global, val_max_global;
        MPI_Allreduce(&val_l2, &val_l2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&val_max, &val_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        val_l2_global = std::sqrt(val_l2_global);

        recordError(val_l2_global, val_max_global, spec.nx * spec.ny * spec.nz);
      }

      HeatPtr heat;
      DiscVectorPtr u_vec;
  };  


    class HeatMMSMultiConvergenceTester : public StandardDiscSetupMulti,
                                          public SlopeCalculator,
                                          public testing::Test
  {
    protected:
      using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

      HeatMMSMultiConvergenceTester()
      {
        setup(3, 1);
      }

      using StandardDiscSetupMulti::setup;

      virtual void setup(const int quad_degree, int sol_degree, const std::vector<Mesh::MeshSpec>& spec,
                         const linear_system::LargeMatrixOptsPetsc& matrix_opts,
                         const std::vector<bool>& is_surf_dirichlet = {true, true, true, true, true, true, true, true, true, true, true})
      {
        StandardDiscSetupMulti::setup(quad_degree, sol_degree, spec, is_surf_dirichlet);

        heat        = std::make_shared<Heat::HeatEquation>(disc);
        u_vec       = makeDiscVector(disc);
        u_solve_vec = makeDiscVector(disc);
        res_vec     = makeDiscVector(disc);

        auto num_dofs     = disc->getDofNumbering()->getNumOwnedDofs();
        this->matrix_opts = matrix_opts;
        auto sparsity     = std::make_shared<linear_system::SparsityPatternMesh>(mesh);
        mat               = std::make_shared<linear_system::LargeMatrixPetsc>(num_dofs, num_dofs, matrix_opts, sparsity);
        assembler         = std::make_shared<linear_system::Assembler>(disc, mat);
      }

      template <typename Tex, typename Tderiv, typename Tsrc>
      void setSolution(Tex ex_sol, Tderiv deriv, Tsrc src, const Heat::VolumeGroupParams& params = Heat::VolumeGroupParams{1, 1, 1})
      {
        for (int i=0; i < disc->getNumVolDiscs(); ++i)
        {
          heat->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
          heat->addVolumeGroupParams(params);
        }

        for (int i=0; i < disc->getNumBCSurfDiscs(); ++i)
        {
          auto surf = disc->getBCSurfDisc(i);
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
        std::vector<DofInt> owned_to_local_dofs;
        mesh->getOwnedLocalDofInfo(owned_to_local_dofs);
        int num_dofs = u_ex->getNumDofs();
        auto err     = makeDiscVector(disc);
        auto tmp     = makeDiscVector(disc);

        auto& u_ex_vec  = u_ex->getVector();
        auto& u_sol_vec = u_sol->getVector();
        auto& err_vec   = err->getVector();
        for (int i=0; i < num_dofs; ++i)
          err_vec[i] = u_ex_vec[i] - u_sol_vec[i];          

        err->markVectorModified();

        heat->applyMassMatrix(err, tmp);
        if (!tmp->isVectorCurrent())
          tmp->syncArrayToVector();
        
        Real val_l2 = 0, val_max = 0;
        auto& tmp_vec = tmp->getVector();
        for (auto i : owned_to_local_dofs)
        {
          val_l2 += err_vec[i] * tmp_vec[i];
          val_max = std::max(val_max, std::abs(err_vec[i]));
        }
        double val_l2_global, val_max_global;
        MPI_Allreduce(&val_l2, &val_l2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&val_max, &val_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        val_l2_global = std::sqrt(val_l2_global);

        //int nelems = specs[0].nx * (specs[0].ny + specs[1].ny) + specs[0].nz;
        Real h = 1.0/specs[0].nx;
        recordError(val_l2_global, val_max_global, h);
      }

      HeatPtr heat;
      DiscVectorPtr u_vec;
      DiscVectorPtr u_solve_vec;
      DiscVectorPtr res_vec;
      linear_system::LargeMatrixOptsPetsc matrix_opts;
      linear_system::LargeMatrixPtr mat;
      linear_system::AssemblerPtr assembler;
  };


  linear_system::LargeMatrixOptsPetsc get_options()
  {
    linear_system::LargeMatrixOptsPetsc opts;
    opts.is_structurally_symmetric = false;
    opts.is_value_symmetric        = false;
    opts.factor_in_place           = false;
    opts.petsc_opts["ksp_atol"] = "1e-12";
    opts.petsc_opts["ksp_rtol"] = "1e-50";
    opts.petsc_opts["ksp_monitor"] = "";

    return opts;
  }
}

TEST(VerificationTest, Foo)
{
  EXPECT_TRUE(true);
}

TEST_F(HeatMMSConvergenceTester, Exponential)
{
  Real kappa = 2;
  Heat::VolumeGroupParams params{kappa, 3, 4};
  for (int sol_degree=1; sol_degree <=2; ++sol_degree)
  {
    auto opts = get_options();
    auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { return std::exp(x + y + z); };

    auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> std::array<Real, 3>
                       { return {kappa * std::exp(x + y + z), kappa * std::exp(x + y + z), kappa * std::exp(x + y + z)}; };
    auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return -3 * kappa * std::exp(x + y + z); };

    int nmeshes = 4;
    int nelem_start = 3;
    for (int i=0; i < nmeshes; ++i)
    {
      int nelem = nelem_start * std::pow(2, i);
      auto meshspec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, nelem, nelem, nelem);

      std::cout << "mesh " << i << " with " << std::pow(nelem, 3) << " elements" << std::endl;
      setup(2*sol_degree, sol_degree, meshspec, opts);
      setSolution(ex_sol_l, deriv_l, src_func_l, params);
      int num_dofs = u_vec->getNumDofs();
      mesh->getFieldDataManager().attachVector(u_vec, "solution");
      mesh->writeVtkFiles(std::string("mesh") + std::to_string(i));

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

    // for some reason, p=2 converges at a rate of 4
    EXPECT_NEAR(slopes[slopes.size()-1][0], 2*sol_degree, 0.1);
    resetErrors();
  }
}


TEST_F(UnsteadyHeatMMSConvergenceTester, CrankNicolsonExponential)
{
  //Real kappa = 2;
  //Heat::VolumeGroupParams params{kappa, 3, 4};

  Real kappa = 1;
  Heat::VolumeGroupParams params{kappa, 1, 1};


  timesolvers::TimeStepperOpts opts;
  opts.t_start = 0.0;
  opts.t_end   = 0.5;
  opts.delta_t = 0.1; 
  opts.mat_type = linear_system::LargeMatrixType::Petsc;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>(get_options());
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 5;  //TODO: test 1

  for (int sol_degree=1; sol_degree <= 1; ++sol_degree)
  {

    auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { 
                          return std::exp(x + y + z + t); 
                          //return t + 1;
                        };

    auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> std::array<Real, 3>
                       { 
                         return {kappa * std::exp(x + y + z + t), kappa * std::exp(x + y + z + t), kappa * std::exp(x + y + z + t)}; 
                         //return {0, 0, 0}; 
                        };
    auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { 
                            return -2 * kappa * std::exp(x + y + z + t);
                          };

    int nmeshes = 1; //3;
    int nelem_start = 81; //3;
    for (int i=0; i < nmeshes; ++i)
    {
      int nelem = nelem_start * std::pow(2, i);
      auto meshspec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, nelem, nelem, nelem);

      std::cout << "mesh " << i << " with " << std::pow(nelem, 3) << " elements" << std::endl;
      setup(2*sol_degree, sol_degree, meshspec, {false, false, false, false, false, false});
      setSolution(ex_sol_l, deriv_l, src_func_l, params);

      mesh->getFieldDataManager().attachVector(u_vec, "solution");
      mesh->writeVtkFiles(std::string("mesh") + std::to_string(i));

      timesolvers::CrankNicolson crank(heat, u_vec, opts);
      crank.solve();

      auto u_ex_vec = makeDiscVector(disc);
      u_ex_vec->setFunc([&](Real x, Real y, Real z){return ex_sol_l(x, y, z, opts.t_end);});
      
      computeErrorNorm(u_ex_vec, u_vec);
      opts.delta_t *= 0.5;
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

    // for some reason, p=2 converges at a rate of 4
    EXPECT_NEAR(slopes[slopes.size()-1][0], 2*sol_degree, 0.1);
    resetErrors();
  }
}


TEST_F(HeatMMSMultiConvergenceTester, Exponential)
{
  Real kappa = 2;
  Heat::VolumeGroupParams params{kappa, 3, 4};
  for (int sol_degree=1; sol_degree <=2; ++sol_degree)
  {
    auto opts = get_options();
    auto ex_sol_l = [&] (Real x, Real y, Real z, Real t) -> Real
                        { return std::exp(x + y + z); };

    auto deriv_l = [&] (Real x, Real y, Real z, Real t) -> std::array<Real, 3>
                       { return {kappa * std::exp(x + y + z), kappa * std::exp(x + y + z), kappa * std::exp(x + y + z)}; };
    auto src_func_l = [&] (Real x, Real y, Real z, Real t) -> Real
                          { return -3 * kappa * std::exp(x + y + z); };

    int nmeshes = 4;
    int nelem_start = 3;
    for (int i=0; i < nmeshes; ++i)
    {
      int nelem = nelem_start * std::pow(2, i);
      auto meshspec1 = Mesh::getMeshSpec(0, 1, 0, 0.5, 0, 1, nelem, nelem, nelem);
      auto meshspec2 = Mesh::getMeshSpec(0, 1, 0.5, 1, 0, 1, nelem, (nelem + 1), nelem);

      std::cout << "mesh " << i << " with " << std::pow(nelem, 3) << " elements" << std::endl;
      setup(2*sol_degree, sol_degree, {meshspec1, meshspec2}, opts);
      setSolution(ex_sol_l, deriv_l, src_func_l, params);
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
    auto h_vals = getHValues();

    for (unsigned int i=0; i < slopes.size(); ++i)
    {
      Real ratio_l2 = 0, ratio_max = 0;
      if (i > 0)
      {
        ratio_l2 = errors_l2[i-1] / errors_l2[i];
        ratio_max = errors_max[i-1] / errors_max[i];
      }

      std::cout << "mesh " << i << ": slope_l2 = " << slopes[i][0] << ", slope_max = " << slopes[i][1] << ", ratio_l2 = " << ratio_l2 << ", ratio_max = " << ratio_max  << ", h_value = " << h_vals[i] << std::endl;
    }

    // for some reason, p=2 converges at a rate of 4
    EXPECT_NEAR(slopes[slopes.size()-1][0], 2*sol_degree, 0.3);
    resetErrors();
  }
}