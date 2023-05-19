#include "file/EpwReader.h"
#include <iostream>

int main(int argc, char* argv[])
{
  if (argc != 2)
    std::cerr << "Usage: " << argv[0] << " fname.epw" << std::endl;

  EPWReader reader(argv[1]);

  return 0;
}
