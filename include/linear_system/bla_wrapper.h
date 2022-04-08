#ifndef BLA_WRAPPER_H
#define BLA_WRAPPER_H

#include "lapacke.h"
#include "cblas.h"
#include <vector>

namespace linear_system {

// Note: all matrices must be column major!

//-----------------------------------------------------------------------------
// GEMV
void gemv(char trans, blasint m, blasint n, double alpha, const double* A, blasint lda, 
          const double* x, blasint incx, double beta, double* y, blasint incy);

void gemv(char trans, blasint m, blasint n, float alpha, const float* A, blasint lda, 
          const float* x, blasint incx, float beta, float* y, blasint incy);

template <typename T>
void gemv(char trans, blasint m, blasint n, T alpha, const std::vector<T>& A, const std::vector<T>& x, T beta, std::vector<T>& y);

extern template
void gemv<double>(char trans, blasint m, blasint n, double alpha, const std::vector<double>& A, const std::vector<double>& x, double beta, std::vector<double>& y);

extern template
void gemv<float>(char trans, blasint m, blasint n, float alpha, const std::vector<float>& A, const std::vector<float>& x, float beta, std::vector<float>& y);


//-----------------------------------------------------------------------------
// SYMV
void symv(char uplo, blasint n, double alpha, const double* A, blasint lda, 
          const double* x, blasint incx, double beta, double* y, blasint incy);

void symv(char uplo, blasint n, float alpha, const float* A, blasint lda, 
          const float* x, blasint incx, float beta, float* y, blasint incy);

template <typename T>
void symv(char uplo, blasint n, T alpha, const std::vector<T>& A, const std::vector<T>& x, T beta, std::vector<T>& y);

extern template
void symv<double>(char uplo, blasint n, double alpha, const std::vector<double>& A, const std::vector<double>& x, double beta, std::vector<double>& y);

extern template
void symv<float>(char uplo, blasint n, float alpha, const std::vector<float>& A, const std::vector<float>& x, float beta, std::vector<float>& y);


//-----------------------------------------------------------------------------
// GESV

lapack_int gesv(lapack_int n, lapack_int nrhs, double* a, lapack_int lda, lapack_int* ipiv, double* b, lapack_int ldb);

lapack_int gesv(lapack_int n, lapack_int nrhs, float* a, lapack_int lda, lapack_int* ipiv, float* b, lapack_int ldb);

template <typename T>
lapack_int gesv(lapack_int n, std::vector<T>& a, std::vector<lapack_int>& ipiv, std::vector<T>& b);

extern template
lapack_int gesv<double>(lapack_int n, std::vector<double>& a, std::vector<lapack_int>& ipiv, std::vector<double>& b);

extern template
lapack_int gesv<float>(lapack_int n, std::vector<float>& a, std::vector<lapack_int>& ipiv, std::vector<float>& b);

//-----------------------------------------------------------------------------
// GETRF

lapack_int getrf(lapack_int m, lapack_int n, double* a, lapack_int lda, lapack_int* ipiv);

lapack_int getrf(lapack_int m, lapack_int n, float* a, lapack_int lda, lapack_int* ipiv);

template <typename T>
lapack_int getrf(lapack_int m, lapack_int n, std::vector<T>& a, std::vector<lapack_int>& ipiv);

extern template
lapack_int getrf<double>(lapack_int m, lapack_int n, std::vector<double>& a, std::vector<lapack_int>& ipiv);

extern template
lapack_int getrf<float>(lapack_int m, lapack_int n, std::vector<float>& a, std::vector<lapack_int>& ipiv);


//-----------------------------------------------------------------------------
// GETRS

lapack_int getrs(char trans, lapack_int n, lapack_int nrhs, double* a, lapack_int lda, const lapack_int* ipiv, double* b, lapack_int ldb);

lapack_int getrs(char trans, lapack_int n, lapack_int nrhs, float* a, lapack_int lda, const lapack_int* ipiv, float* b, lapack_int ldb);

template <typename T>
lapack_int getrs(char trans, lapack_int n, std::vector<T>& a, const std::vector<lapack_int>& ipiv, std::vector<T>& b);

extern template
lapack_int getrs<double>(char trans, lapack_int n, std::vector<double>& a, const std::vector<lapack_int>& ipiv, std::vector<double>& b);

extern template
lapack_int getrs<float>(char trans, lapack_int n, std::vector<float>& a, const std::vector<lapack_int>& ipiv, std::vector<float>& b);


//-----------------------------------------------------------------------------
// POTRF

lapack_int potrf(char uplo, lapack_int n, double* a, lapack_int lda);

lapack_int potrf(char uplo, lapack_int n, float* a, lapack_int lda);

template <typename T>
lapack_int potrf(char uplo, lapack_int n, std::vector<T>& a);

extern template
lapack_int potrf<double>(char uplo, lapack_int n, std::vector<double>& a);

extern template
lapack_int potrf<float>(char uplo, lapack_int n, std::vector<float>& a);

//-----------------------------------------------------------------------------
// POTRS

lapack_int potrs(char uplo, lapack_int n, lapack_int nrhs, const double* a, lapack_int lda, double* b, lapack_int ldb);

lapack_int potrs(char uplo, lapack_int n, lapack_int nrhs, const float* a, lapack_int lda, float* b, lapack_int ldb);

template <typename T>
lapack_int potrs(char uplo, lapack_int n, const std::vector<T>& a, std::vector<T>& b);

extern template
lapack_int potrs<double>(char uplo, lapack_int n, const std::vector<double>& a, std::vector<double>& b);

extern template
lapack_int potrs<float>(char uplo, lapack_int n, const std::vector<float>& a, std::vector<float>& b);

}  // namespace

#endif