#ifndef UTIL_MEMORY_H
#define UTIL_MEMORY_H

#include <memory>

template <typename T, typename Deleter>
std::shared_ptr<T> makeSharedWithDeleter(T* p, Deleter& deleter)
{
  return std::shared_ptr<T>(p, deleter);
}

#endif