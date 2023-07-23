#include "discretization/NeumannBC_defs.h"
#include "discretization/DirichletBC_defs.h"
#include "discretization/discretization.h"
#include "discretization/surface_discretization.h"
#include "discretization/volume_discretization.h"
#include "error_handling.h"
#include "initialization.h"
#include "linear_system/large_matrix_petsc.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/heat/HeatEquation.h"
#include "time_solver/crank_nicolson.h"
#include "time_solver/time_stepper_opts.h"
#include "time_solver/timestep_controller.h"
#include "utils/string_utils.h"
#include "mesh/mesh_generator.h"
#include "mesh/mesh.h"

struct Opts
{
  int numel_x = 0;
  int numel_y = 0;
  int numel_z = 0;
  double delta_t = -1;
  double tmax    = 0;
};

Opts parseArguments(int argc, char* argv[])
{
  Opts opts;
  Parser parser;

  if (argc != 6)
  {
    std::cout << "Usage: " << argv[0] << "numel_x numel_y numel_z delta_t tmax" << std::endl;
    throw std::runtime_error("incorrect number of command line arguments");
  }

  opts.numel_x = parser.get<int>(argv[1]);
  opts.numel_y = parser.get<int>(argv[2]);
  opts.numel_z = parser.get<int>(argv[3]);
  opts.delta_t = parser.get<double>(argv[4]);
  opts.tmax    = parser.get<double>(argv[5]);

  return opts;
}

std::shared_ptr<Mesh::MeshCG> createMesh(const Opts& opts)
{
  int sol_degree = 1;
  std::vector<bool> is_surf_dirichlet = {true, false, false, false, false};

  Mesh::MeshSpec meshspec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, opts.numel_x, opts.numel_y, opts.numel_z);
  
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  auto m = Mesh::make_parallel_mesh(meshspec, comm_size, &(Mesh::identity));

  // make volume groups
  Mesh::MeshEntityGroupSpec volume_group("volume0");
  volume_group.addModelEntity(Mesh::ModelEntitySpec(3, 0));
  std::vector<Mesh::MeshEntityGroupSpec> volume_groups{volume_group};

  // make surface groups
  std::vector<Mesh::MeshEntityGroupSpec> surface_groups;
  for (int i=0; i < 6; ++i)
  {
    surface_groups.emplace_back(std::string("surface") + std::to_string(i));
    surface_groups.back().addModelEntity(Mesh::ModelEntitySpec(2, i), Mesh::ModelEntitySpec(3, 0));
    surface_groups.back().setIsDirichlet(is_surf_dirichlet[i]);
  }

  std::vector<Mesh::MeshEntityGroupSpec> other_surfaces;
  auto mesh = Mesh::createMeshCG(m, volume_groups, surface_groups,
                                 other_surfaces, sol_degree, 1);

  return mesh;
}

std::shared_ptr<Heat::HeatEquation> createPhysics(std::shared_ptr<Mesh::MeshCG> mesh)
{
  int quad_degree = 3;
  auto disc = std::make_shared<Discretization>(mesh, quad_degree, quad_degree);
  auto heat = std::make_shared<Heat::HeatEquation>(disc);

  return heat;
}

void setSolution(std::shared_ptr<Heat::HeatEquation> heat, DiscVectorPtr sol)
{
  Real rho = 64, Cp = 2020, kappa = 0.039;
  double pi = std::atan2(0, -1);

  Heat::VolumeGroupParams params = {kappa, rho, Cp};
  auto ex_sol = [=](const Real& x, const Real& y, const Real& z, const Real& t)
  {
    return x*x + y*y + z*z + std::sin(pi * t/86400);
  };

  auto ex_sol_dt = [=](const Real& x, const Real& y, const Real& z, const Real& t)
  {
    return pi * std::cos(pi * t / 86400) / 86400;
  };

  auto deriv = [=](const Real& x, const Real& y, const Real& z, const Real& t)
  {
    return std::array<Real, 3>{2*x*kappa, 2*y*kappa, 2*z*kappa};
  };

  auto src = [=](const Real& x, const Real& y, const Real& z, const Real& t)
  {
    return rho * Cp * pi * std::cos(pi * t / 86400) / 86400 - 8*kappa;
  };


  auto disc = heat->getDiscretization();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    heat->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
    heat->addVolumeGroupParams(params);
  }

  for (int i=0; i < disc->getNumSurfDiscs(); ++i)
  {
    auto surf = disc->getSurfDisc(i);
    if (surf->getIsDirichlet())
      heat->addDirichletBC(makeDirichletBCMMS(surf, ex_sol, ex_sol_dt));
    else
      heat->addNeumannBC(makeNeumannBCMMS(surf, deriv));
  }

  auto f = [&](Real x, Real y, Real z)
            { return ex_sol(x, y, z, 0); };
  sol->setFunc(f);
}

timesolvers::TimeStepperOpts getTimeStepperOpts(const Opts& input_opts)
{
  timesolvers::TimeStepperOpts opts;

  opts.t_start = 0;
  opts.t_end   = input_opts.tmax; 
  opts.nonlinear_abs_tol = 1e-9; 
  opts.nonlinear_rel_tol = 1e-9; 
  opts.nonlinear_itermax = 30; 
  opts.vis_output_freq   = -1; 

  opts.timestep_controller = std::make_shared<timesolvers::TimestepControllerConstant>(input_opts.delta_t);
  opts.mat_type = linear_system::LargeMatrixType::Petsc;


  auto matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>();
  matrix_opts->is_structurally_symmetric = true;
  opts.solve_auxiliary_equations_combined_system = true;
  opts.precompute_linear_jacobian = true;

  matrix_opts->petsc_opts["ksp_atol"] = "1e-12"; //input_vals.at("linear_abs_tol");
  matrix_opts->petsc_opts["ksp_rtol"] = "1e-12"; //input_vals.at("linear_rel_tol");
  matrix_opts->petsc_opts["ksp_monitor"] = "";

  if (commSize(MPI_COMM_WORLD) > 1)
  {
    matrix_opts->petsc_opts["ksp_type"] = "cg";
    //matrix_opts->petsc_opts["pc_type"] = "jacobi";

    matrix_opts->petsc_opts["pc_type"] = "asm";
    matrix_opts->petsc_opts["pc_asm_overlap"] = "1";
  }

  opts.matrix_opts = matrix_opts;

  return opts;
}

int getNumDofs(std::shared_ptr<Mesh::MeshCG> mesh)
{
  int num_owned_dofs = mesh->getNumOwnedDofs();
  int num_total_dofs = 0;
  MPI_Allreduce(&num_owned_dofs, &num_total_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  return num_total_dofs;
}



int main(int argc, char* argv[])
{
  PetscOptionsSetValue(NULL, "-on_error_abort", "");
  //linear_system::setPetscGlobalOption("log_view", "");

  initialize(argc, argv);  

  bool am_i_root = commRank(MPI_COMM_WORLD) == 0;

  // Runs a simple MMS case, solving the heat equation with solution T = x^2 + y^2 + z^2 + t^2
  Opts opts = parseArguments(argc, argv);

  auto mesh = createMesh(opts);
  auto heat = createPhysics(mesh);
  auto sol = makeDiscVector(heat->getDiscretization());
  setSolution(heat, sol);
  auto timestepper_opts = getTimeStepperOpts(opts);
  auto u_aux = makeAuxiliaryEquationsStorage(heat->getAuxEquations());

  int num_total_dofs = getNumDofs(mesh);
  if (am_i_root)
  {
    
    std::cout << "number of dofs = " << num_total_dofs << std::endl;
  }

  timesolvers::CrankNicolson timesolver(heat, sol, u_aux, timestepper_opts);

  // run solver
  double t_start_run = MPI_Wtime();
  mesh->getFieldDataManager().attachVector(sol, "solution");
  mesh->writeVtkFiles("solution_initial");
  timesolver.solve();
  mesh->writeVtkFiles("solution_final");

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_run = MPI_Wtime();

  if (commRank(MPI_COMM_WORLD) == 0)
  {
    std::cout << "\n\nFinished simple heat run" << std::endl;
    double t_run_elapsed = t_end_run - t_start_run;
    double t_simulated_elapsed = timestepper_opts.t_end - timestepper_opts.t_start;
    std::cout << "simulation took " << t_run_elapsed << " seconds to simulate " << t_simulated_elapsed << " seconds"
              << ", which is " << t_simulated_elapsed/t_run_elapsed << "x realtime" << std::endl;
  }

  return 0;
}