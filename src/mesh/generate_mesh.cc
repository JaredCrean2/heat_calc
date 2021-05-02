#include <iostream>
#include <stdio.h>
#include "mesh/mesh_generator.h"
#include "mpi.h"
#include "PCU.h"

#include "lionPrint.h"

int main(int argc, char* argv[])
{
  lion_set_verbosity(10);

  std::cout << "argc = " << argc << std::endl;
  if (argc != 4 && argc != 10)
  {
    std::cerr << "Usage: " << argv[0] << " nx ny nz [xmin] [xmax] [ymin] [ymax] [zmin] [zmax]" << std::endl;
    return 1;
  }

  MPI_Init(&argc, &argv);
  PCU_Comm_Init();


  int nx, ny, nz;
  double xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1;
  sscanf(argv[1], "%d", &nx);
  sscanf(argv[2], "%d", &ny);
  sscanf(argv[3], "%d", &nz);
  if (argc == 10)
  {
    sscanf(argv[4], "%lf", &xmin);
    sscanf(argv[5], "%lf", &xmax);
    sscanf(argv[6], "%lf", &ymin);
    sscanf(argv[7], "%lf", &ymax);
    sscanf(argv[8], "%lf", &zmin);
    sscanf(argv[9], "%lf", &zmax);
  }

  Mesh::MeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.xmin = xmin;
  spec.xmax = xmax;
  spec.ymin = ymin;
  spec.ymax = ymax;
  spec.zmin = zmin;
  spec.zmax = zmax;

  auto generator = make_mesh_generator(spec, &(Mesh::identity));
  apf::Mesh2* mesh = generator.generate();

  mesh->writeNative("./mesh_.smb");

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
  

  return 0;
}
