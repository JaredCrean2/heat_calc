#include <vector>
#include <string>
#include <iostream>

using ValType = size_t;
const int BlockSize = 512;

size_t parseSize(const std::string& size, const std::string& suffix)
{
  size_t multiplier;
  if (suffix == "B")
  {
    multiplier = 1;
  } else if (suffix == "K")
  {
    multiplier = 1024;
  } else if (suffix == "M")
  {
    multiplier = 1024*1024;
  } else if (suffix == "G")
  {
    multiplier = 1024 * 1024 * 1024;
  } else
    throw std::runtime_error(std::string("unrecognized suffix ") + suffix);

  
  size_t size_val = std::stoi(size);

  return size_val * multiplier;
}

void allocateMem(size_t size_bytes)
{
  size_t block_size_bytes = sizeof(ValType) * BlockSize;

  int nblocks_whole = size_bytes / block_size_bytes;
  size_t rem_size_bytes = size_bytes - nblocks_whole * block_size_bytes;
  size_t rem_size = rem_size_bytes / sizeof(ValType);

  int nblocks_rem = rem_size > 0;

  std::vector<std::vector<ValType>> data(nblocks_whole + nblocks_rem);

  for (int i=0; i < nblocks_whole; ++i)
  {
    std::cout << "allocating block " << i << " / " << nblocks_whole << std::endl;

    data[i].resize(BlockSize);

    // Some operating system allocate on first touch, so write data
    // to the vector to make sure it really gets allocated
    for (auto& val : data[i])
      val = 1;
  }

  if (nblocks_rem > 0)
  {
    std::cout << "allocating remainder block of size " << rem_size_bytes << " bytes" << std::endl;
    data[nblocks_whole].resize(rem_size);

    for (auto& val : data[nblocks_whole])
      val = 1;
  }
}

size_t parseArgs(int argc, char* argv[])
{
  if (argc != 2 && argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " size [suffix]" << std::endl;
    throw std::runtime_error("unable to parse args");
  }

  std::string size = argv[1];
  std::string suffix = "B";
  if (argc == 3)
    suffix = argv[2];

  return parseSize(size, suffix);
}

int main(int argc, char* argv[])
{

  size_t size = parseArgs(argc, argv);
  allocateMem(size);

  return 0;
}