#include "utils/initialization.h"
#include "petscksp.h"

void initialize(int& argc, ArgvType argv)
{
  PetscInitialize(&argc, &(argv), nullptr, nullptr);
}