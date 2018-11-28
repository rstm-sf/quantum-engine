// Copyright (C) 2018 Rustam Sayfutdinov (rstm.sf@gmail.com)
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef QENGINE_UTILS_MATRIX_H_
#define QENGINE_UTILS_MATRIX_H_

#include <complex>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include <gsl/gsl_assert>

namespace qengine {
inline namespace util {

template <typename T>
class Matrix {
public:
  Matrix<T>();
  virtual ~Matrix<T>();
  Matrix<T>(const Matrix<T>&);
  Matrix<T>(Matrix<T>&&);
  Matrix<T>& operator=(const Matrix<T>&);
  Matrix<T>& operator=(Matrix<T>&&);

  Matrix<T>(uint64_t ncols, uint64_t nrows, const std::vector<T>& vals);
  Matrix<T>(
      uint64_t ncols, uint64_t nrows, const std::initializer_list<T>& vals);
  Matrix<T>(uint64_t ncols, uint64_t nrows);

  T& operator()(uint64_t i, uint64_t j);
  T operator()(uint64_t i, uint64_t j) const;
  Matrix<T> operator+(const Matrix<T>& A);
  Matrix<T> operator-(const Matrix<T>& A);
  Matrix<T> operator+(const Matrix<T>& A) const;
  Matrix<T> operator-(const Matrix<T>& A) const;
  Matrix<T>& operator+=(const Matrix<T>& A);
  Matrix<T>& operator-=(const Matrix<T>& A);
  Matrix<T>& operator*=(const Matrix<T>& A);

  T trace() const;
  Matrix<T> transpose() const;
  Matrix<T> dagger() const;
  Matrix<T> conjugate() const;
  Matrix<T> tensor_times(const Matrix<T> & A) const;

  uint64_t ncols() const;
  uint64_t nrows() const;
  std::vector<T> vals() const;
  uint64_t size() const;

  template <typename T1>
  friend std::ostream& operator<<(std::ostream& out, const Matrix<T1>& A);

  template <typename T1>
  friend std::istream& operator>>(std::istream& in, const Matrix<T1>& A);

  template <typename T1, typename T2>
  friend Matrix<T1> operator*(T2 alpha, const Matrix<T1>& A);

  template <typename T1, typename T2>
  friend Matrix<T1> operator*(const Matrix<T1>& A, T2 alpha);

  template <typename T1, typename T2>
  friend std::vector<T2> operator*(
      const Matrix<T1>& A, const std::vector<T2>& b);

  template <typename T1, typename T2>
  friend std::vector<T1> operator*(
      const std::vector<T1>& b, const Matrix<T2>& A);

  template <typename T1>
  friend Matrix<T1> operator*(const Matrix<T1>& A, const Matrix<T1>& B);

  template <typename T1>
  friend bool operator==(const Matrix<T1>& A, const Matrix<T1>& B);

  template <typename T1>
  friend bool operator!=(const Matrix<T1>& A, const Matrix<T1>& B);

protected:
  uint64_t ncols_;
  uint64_t nrows_;
  std::vector<T> vals_;
};

template <typename T>
Matrix<T>::Matrix() = default;

template <typename T>
Matrix<T>::~Matrix() = default;

template <typename T>
Matrix<T>::Matrix(const Matrix<T>&) = default;

template <typename T>
Matrix<T>::Matrix(Matrix<T>&&) = default;

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>&) = default;

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&&) = default;

template <typename T>
Matrix<T>::Matrix(uint64_t nrows, uint64_t ncols, const std::vector<T>& vals)
  : nrows_{nrows}, ncols_{ncols}, vals_{vals} {
  Expects(nrows_ * ncols_ == this->size());
}

template <typename T>
Matrix<T>::Matrix(
    uint64_t nrows, uint64_t ncols, const std::initializer_list<T>& vals)
  : nrows_{nrows}, ncols_{ncols}, vals_{vals} {
  Expects(nrows_ * ncols_ == this->size());
}

template <typename T>
Matrix<T>::Matrix(uint64_t nrows, uint64_t ncols)
  : nrows_{nrows}, ncols_{ncols}, vals_(nrows * ncols) {}

template <typename T>
T& Matrix<T>::operator()(uint64_t i, uint64_t j) {
  return vals_[i + j * nrows_];
}

template <typename T>
T Matrix<T>::operator()(uint64_t i, uint64_t j) const {
  return vals_[i + j * nrows_];
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) {
  Expects(ncols_ == A.ncols_ && nrows_ == A.nrows_);

  Matrix<T> B(*this);
  for (uint64_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] += A.vals_[i];
  return B;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) {
  Expects(ncols_ == A.ncols_ && nrows_ == A.nrows_);

  Matrix<T> B(*this);
  for (uint64_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] -= A.vals_[i];
  return B;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) const {
  Expects(ncols_ == A.ncols_ && nrows_ == A.nrows_);

  Matrix<T> B(*this);
  for (uint64_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] += A.vals_[i];
  return B;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) const {
  Expects(ncols_ == A.ncols_ && nrows_ == A.nrows_);

  Matrix<T> B(*this);
  for (uint64_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] -= A.vals_[i];
  return B;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& A) {
  Expects(ncols_ == A.ncols_ && nrows_ == A.nrows_);

  for (uint64_t i = 0; i < vals_.size(); ++i)
    vals_[i] += A.vals_[i];
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& A) {
  Expects(ncols_ == A.ncols_ && nrows_ == A.nrows_);

  for (uint64_t i = 0; i < vals_.size(); ++i)
    vals_[i] -= A.vals_[i];
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& A) { return *this * A; }

template <class T>
T Matrix<T>::trace() const {
  Expects(nrows_ == ncols_);

  T sum{};
  for (uint64_t i = 0; i < nrows_; ++i)
    sum += vals_[i + i * nrows_];
  return sum;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> temp(ncols_, nrows_);
  for (uint64_t i = 0; i < ncols_; ++i)
    for (uint64_t j = 0; j < nrows_; ++j)
      temp.vals_[i + j * ncols_] = vals_[j + i * nrows_];
  return temp;
}

template <typename T>
Matrix<T> Matrix<T>::dagger() const {
  Matrix<T> temp(ncols_, nrows_);
  for (uint64_t i = 0; i < ncols_; ++i)
    for (uint64_t j = 0; j < nrows_; ++j)
      temp.vals_[i + j * ncols_] = std::conj(vals_[j + i * nrows_]);
  return temp;
}

template <typename T>
Matrix<T> Matrix<T>::conjugate() const {
  Matrix<T> temp(nrows_, ncols_);
  for (uint64_t i = 0; i < nrows_; ++i)
    for (uint64_t j = 0; j < ncols_; ++j)
      temp.vals_[i + j * nrows_] = std::conj(vals_[i + j * nrows_]);
  return temp;
}

template <class T>
Matrix<T> Matrix<T>::tensor_times(const Matrix<T> & A) const {
  Matrix<T> C(nrows_ * A.nrows_, ncols_ * A.ncols_);
  for (uint64_t q = 0; q < A.ncols_; ++q)
    for (uint64_t j = 0; j < ncols_; ++j)
      for (uint64_t i = 0; i < nrows_; ++i)
        for (uint64_t p = 0; p < A.nrows_; ++p) {
          uint64_t n = i * A.nrows_ + p;
          uint64_t m = j * A.ncols_ + q;
          C(n, m) = vals_[i + j * nrows_] * A(p, q);
        }
  return C;
}

template <typename T>
uint64_t Matrix<T>::ncols() const { return ncols_; }

template <typename T>
uint64_t Matrix<T>::nrows() const { return nrows_; }

template <typename T>
std::vector<T> Matrix<T>::vals() const { return vals_; }

template <typename T>
uint64_t Matrix<T>::size() const {
  return static_cast<uint64_t>(vals_.size());
}

template <typename T1>
std::ostream& operator<<(std::ostream& out, const Matrix<T1>& A) {
  for (uint64_t i = 0; i < A.nrows_; ++i) {
    for (uint64_t j = 0; j < A.ncols_; ++j)
      out << A.vals_[i + j * A.nrows_] << "\t";
    out << "\n";
  }
  return out;
}

template <typename T>
std::istream& operator>>(std::istream& in, const Matrix<T>& A) {
  for (uint64_t j = 0; j < A.ncols_; ++j)
    for (uint64_t i = 0; i < A.nrows_; ++i)
      in >> A.vals_[i + j * A.nrows_];
  return in;
}

template <typename T1, typename T2>
Matrix<T1> operator*(T2 alpha, const Matrix<T1>& A) {
  Matrix<T1> B(A);
  T1 beta = static_cast<T1>(alpha);
  for (auto & v : B.vals_)
    v *= beta;
  return B;
}

template <typename T1, typename T2>
Matrix<T1> operator*(const Matrix<T1>& A, T2 alpha) { return alpha * A; }

template <typename T1, typename T2>
std::vector<T2> operator*(const Matrix<T1>& A, const std::vector<T2>& b) {
  Expects(A.ncols_ == b.size());

  std::vector<T2> res(A.nrows_);
  for (uint64_t j = 0; j < A.ncols_; ++j)
    for (uint64_t i = 0; i < A.nrows_; ++i)
      res[i] += A.vals_[i + j * A.nrows_] * b[j];
  return res;
}

template <typename T1, typename T2>
std::vector<T1> operator*(const std::vector<T1>& b, const Matrix<T2>& A) {
  Expects(A.nrows_ == b.size());

  std::vector<T1> res(A.ncols_);
  for (uint64_t j = 0; j < A.ncols_; ++j)
    for (uint64_t i = 0; i < A.nrows_; ++i)
      res[j] += b[i] * A.vals_[i + j * A.nrows_];
  return res;
}

template <typename T1>
Matrix<T1> operator*(const Matrix<T1>& A, const Matrix<T1>& B) {
  Expects(A.ncols_ == B.nrows_);

  Matrix<T1> C(A.nrows_, B.ncols_);
  // column-major order: A(i, j) = A.val[i + j * nrows_]
  for(uint64_t j = 0; j < B.ncols_; ++j)
    for(uint64_t k = 0; k < A.ncols_; ++k)
      for(uint64_t i = 0; i < A.nrows_; ++i)
        C(i, j) += A(i, k) * B(k, j);

  return C;
}

template <typename T1>
bool operator==(const Matrix<T1>& A, const Matrix<T1>& B) {
  if (A.ncols_ != B.ncols_ || A.nrows_ != B.nrows_ )
    return false;

  return A.vals_ == B.vals_;
}

template <typename T1>
bool operator!=(const Matrix<T1>& A, const Matrix<T1>& B) {
  return !(A.vals_ == B.vals_);
}

} // namespace util
} // namespace qengine

#endif // QENGINE_UTILS_MATRIX_H_