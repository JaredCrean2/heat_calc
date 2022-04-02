#include "linear_system/bla_wrapper.h"

#include <string>
#include <stdexcept>
#include <cassert>
#include <iostream>

namespace linear_system {

namespace {

void check_error(const char* fname, int ret_val)
{
  if (ret_val != 0)
    throw std::runtime_error(std::string("Function ") + fname + " return value " + std::to_string(ret_val));
}

template <typename T>
lapack_int size(const std::vector<T>& v)
{
  return v.size();
}

}



//-----------------------------------------------------------------------------
// GESV

lapack_int gesv(lapack_int n, lapack_int nrhs, double* a, lapack_int lda, lapack_int* ipiv, double* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_dgesv_work(LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("dgesv", ret_val);

  return ret_val;
}

lapack_int gesv(lapack_int n, lapack_int nrhs, float* a, lapack_int lda, lapack_int* ipiv, float* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_sgesv_work(LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("sgesv", ret_val);

  return ret_val;
}

template <typename T>
lapack_int gesv(lapack_int n, std::vector<T>& a, std::vector<lapack_int>& ipiv, std::vector<T>& b)
{
  assert(n > 0);
  assert(n * n == size(a));
  assert(size(ipiv) == n);
  assert(size(b) % n == 0);
  lapack_int lda = n;
  lapack_int ldb = n;
  lapack_int nrhs = size(b) / n;

  return gesv(n, nrhs, a.data(), lda, ipiv.data(), b.data(), ldb);
}

template
lapack_int gesv<double>(lapack_int n, std::vector<double>& a, std::vector<lapack_int>& ipiv, std::vector<double>& b);

template
lapack_int gesv<float>(lapack_int n, std::vector<float>& a, std::vector<lapack_int>& ipiv, std::vector<float>& b);

//-----------------------------------------------------------------------------
// GETRF

lapack_int getrf(lapack_int m, lapack_int n, double* a, lapack_int lda, lapack_int* ipiv)
{
  auto ret_val = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, m, n, a, lda, ipiv);
  check_error("dgetrf", ret_val);

  return ret_val;
}

lapack_int getrf(lapack_int m, lapack_int n, float* a, lapack_int lda, lapack_int* ipiv)
{
  auto ret_val = LAPACKE_sgetrf_work(LAPACK_COL_MAJOR, m, n, a, lda, ipiv);
  check_error("sgetrf", ret_val);

  return ret_val;
}

template <typename T>
lapack_int getrf(lapack_int m, lapack_int n, std::vector<T>& a, std::vector<lapack_int>& ipiv)
{
  assert(m > 0);
  assert(n > 0);
  assert(m * n == size(a));
  assert(size(ipiv) == std::min(m, n));
  lapack_int lda = m;

  return getrf(m, n, a.data(), lda, ipiv.data());
}

template
lapack_int getrf<double>(lapack_int m, lapack_int n, std::vector<double>& a, std::vector<lapack_int>& ipiv);

template
lapack_int getrf<float>(lapack_int m, lapack_int n, std::vector<float>& a, std::vector<lapack_int>& ipiv);


//-----------------------------------------------------------------------------
// GETRS

lapack_int getrs(char trans, lapack_int n, lapack_int nrhs, double* a, lapack_int lda, const lapack_int* ipiv, double* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_dgetrs_work(LAPACK_COL_MAJOR, trans, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("dgetrs", ret_val);

  return ret_val;
}

lapack_int getrs(char trans, lapack_int n, lapack_int nrhs, float* a, lapack_int lda, const lapack_int* ipiv, float* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_sgetrs_work(LAPACK_COL_MAJOR, trans, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("sgetrs", ret_val);

  return ret_val;
}

template <typename T>
lapack_int getrs(char trans, lapack_int n, std::vector<T>& a, const std::vector<lapack_int>& ipiv, std::vector<T>& b)
{
  assert(n > 0);
  assert(n * n == size(a));
  assert(size(ipiv) == n);
  assert(size(b) % n == 0);
  lapack_int lda = n;
  lapack_int ldb = n;
  lapack_int nrhs = size(b) / n;

  return getrs(trans, n, nrhs, a.data(), lda, ipiv.data(), b.data(), ldb);
}

template
lapack_int getrs<double>(char trans, lapack_int n, std::vector<double>& a, const std::vector<lapack_int>& ipiv, std::vector<double>& b);

template
lapack_int getrs<float>(char trans, lapack_int n, std::vector<float>& a, const std::vector<lapack_int>& ipiv, std::vector<float>& b);


//-----------------------------------------------------------------------------
// POTRF

lapack_int potrf(char uplo, lapack_int n, double* a, lapack_int lda)
{
  auto ret_val = LAPACKE_dpotrf_work(LAPACK_COL_MAJOR, uplo, n, a, lda);
  check_error("dpotrf", ret_val);

  return ret_val;
}

lapack_int potrf(char uplo, lapack_int n, float* a, lapack_int lda)
{
  auto ret_val = LAPACKE_spotrf_work(LAPACK_COL_MAJOR, uplo, n, a, lda);
  check_error("spotrf", ret_val);

  return ret_val;
}

template <typename T>
lapack_int potrf(char uplo, lapack_int n, std::vector<T>& a)
{
  assert(n > 0);
  assert(size(a) == n*n);
  lapack_int lda = n;

  return potrf(uplo, n, a.data(), lda);
}

template
lapack_int potrf<double>(char uplo, lapack_int n, std::vector<double>& a);

template
lapack_int potrf<float>(char uplo, lapack_int n, std::vector<float>& a);


//-----------------------------------------------------------------------------
// POTRS

lapack_int potrs(char uplo, lapack_int n, lapack_int nrhs, const double* a, lapack_int lda, double* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_dpotrs_work(LAPACK_COL_MAJOR, uplo, n, nrhs, a, lda, b, ldb);
  check_error("dpotrs", ret_val);

  return ret_val;
}

lapack_int potrs(char uplo, lapack_int n, lapack_int nrhs, const float* a, lapack_int lda, float* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_spotrs_work(LAPACK_COL_MAJOR, uplo, n, nrhs, a, lda, b, ldb);
  check_error("spotrs", ret_val);

  return ret_val;
}

template <typename T>
lapack_int potrs(char uplo, lapack_int n, const std::vector<T>& a, std::vector<T>& b)
{
  assert( n > 0);
  assert(size(a) == n*n);
  assert(size(b) % n == 0);
  lapack_int lda = n;
  lapack_int ldb = n;
  lapack_int nrhs = size(b) / n;


  return potrs(uplo, n, nrhs, a.data(), lda, b.data(), ldb);
}

template
lapack_int potrs<double>(char uplo, lapack_int n, const std::vector<double>& a, std::vector<double>& b);

template
lapack_int potrs<float>(char uplo, lapack_int n, const std::vector<float>& a, std::vector<float>& b);


} // namespace