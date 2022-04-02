#ifndef BLA_WRAPPER_H
#define BLA_WRAPPER_H

#include "lapacke.h"
#include <vector>

namespace linear_system {

// Note: all matrices must be column major!

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