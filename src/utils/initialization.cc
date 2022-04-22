#include "utils/initialization.h"
#include "petscksp.h"
#include <petscsys.h>

void initialize(int& argc, ArgvType argv)
{
  PetscInitialize(&argc, &(argv), nullptr, nullptr);
}

void finalize()
{
  PetscFinalize();
}