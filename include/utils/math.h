#ifndef UTILS_MATH_H
#define UTILS_MATH_H

#include <array>
#include <ostream>

namespace {
  template <typename T>
  using Vec3 = std::array<T, 3>;
}

template <typename T,  size_t N>
std::array<T, N> operator+(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] + b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
std::array<T, N> operator+(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] + b;

  return c;
}

template <typename T,  size_t N>
std::array<T, N> operator-(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] - b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
std::array<T, N> operator-(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] - b;

  return c;
}

template <typename T,  size_t N>
std::array<T, N> operator*(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] * b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
std::array<T, N> operator*(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] * b;

  return c;
}

template <typename T,  size_t N>
std::array<T, N> operator/(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] / b[i];

  return c;
}


template <typename T,  size_t N, typename T2>
std::array<T, N> operator/(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] / b;

  return c;
}

template <typename T,  size_t N>
T dot(const std::array<T, N>& a, const std::array<T, N>& b)
{
  T val = 0;
  for (int i=0; i < N; ++i)
    val += a[i] * b[i];

  return val;
}

template <typename T,  size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& a)
{
  for (size_t i=0; i < N; ++i)
  {
    os << a[i];
    if (i < N - 1)
      os << ", ";
  }

  return os;
}

template <typename T>
Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
  T c1 =   a[1] * b[2] - a[2] * b[1];
  T c2 = -(a[0] * b[2] - a[2] * b[0]);
  T c3 =   a[0] * b[1] - a[1] * b[0];


  return {c1, c2, c3};
}


#endif