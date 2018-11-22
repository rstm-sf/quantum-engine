#ifndef QENGINE_UTILS_MATRIX_H_
#define QENGINE_UTILS_MATRIX_H_

#include <initializer_list>
#include <vector>

namespace qengine {
inline namespace util {

template <typename T>
class Matrix {
public:
  Matrix<T>();
  ~Matrix<T>();
  Matrix<T>(const Matrix<T>&);
  Matrix<T>(Matrix<T>&&);
  Matrix<T>& operator=(const Matrix<T>&);
  Matrix<T>& operator=(Matrix<T>&&);

  Matrix<T>(size_t ncols, size_t nrows, const std::vector<T>& vals);
  Matrix<T>(size_t ncols, size_t nrows, const std::initializer_list<T>& vals);

  T& operator()(size_t i, size_t j);
  T operator()(size_t i, size_t j) const;
  Matrix<T> operator+(const Matrix<T>& A);
  Matrix<T> operator-(const Matrix<T>& A);
  Matrix<T> operator+(const Matrix<T>& A) const;
  Matrix<T> operator-(const Matrix<T>& A) const;
  Matrix<T>& operator+=(const Matrix<T>& A);
  Matrix<T>& operator-=(const Matrix<T>& A);

  size_t get_ncols() const;
  size_t get_nrows() const;
  std::vector<T> get_vals() const;

  template <typename T1>
  friend std::ostream& operator<<(std::ostream& out, const Matrix<T1>& A);

  template <typename T1>
  friend std::istream& operator>>(std::istream& in, const Matrix<T1>& A);

  template <typename T1, typename T2>
  friend Matrix<T1> operator*(T2 alpha, const Matrix<T1>& A);

  template <typename T1, typename T2>
  friend Matrix<T1> operator*(const Matrix<T1>& A, T2 alpha);

  template <typename T1>
  friend Matrix<T1> operator*(const Matrix<T1>& A, const Matrix<T1>& B);

  template <typename T1>
  friend bool operator==(const Matrix<T1>& A, const Matrix<T1>& B);

  template <typename T1>
  friend bool operator!=(const Matrix<T1>& A, const Matrix<T1>& B);

private:
  size_t ncols_;
  size_t nrows_;
  std::vector<T> vals_;

  Matrix<T>(size_t ncols, size_t nrows);
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
Matrix<T>::Matrix(size_t nrows, size_t ncols, const std::vector<T>& vals)
    : nrows_{nrows}, ncols_{ncols}, vals_{vals} {}

template <typename T>
Matrix<T>::Matrix(
    size_t nrows, size_t ncols, const std::initializer_list<T>& vals)
    : nrows_{nrows}, ncols_{ncols}, vals_{vals} {}

template <typename T>
Matrix<T>::Matrix(size_t nrows, size_t ncols)
    : nrows_{nrows}, ncols_{ncols}, vals_(nrows * ncols) {}

template <typename T>
T& Matrix<T>::operator()(size_t i, size_t j) {
  return vals_[i + j * nrows_];
}

template <typename T>
T Matrix<T>::operator()(size_t i, size_t j) const {
  return vals_[i + j * nrows_];
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) {
  if (ncols_ != A.get_ncols() || nrows_ != A.get_nrows())
    throw std::runtime_error{"error: mult_mat: Matrices are not consistent"};

  Matrix<T> B(*this);
  for (size_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] += A.vals_[i];
  return B;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) {
  if (ncols_ != A.get_ncols() || nrows_ != A.get_nrows())
    throw std::runtime_error{"error: mult_mat: Matrices are not consistent"};

  Matrix<T> B(*this);
  for (size_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] -= A.vals_[i];
  return B;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) const {
  if (ncols_ != A.get_ncols() || nrows_ != A.get_nrows())
    throw std::runtime_error{"error: mult_mat: Matrices are not consistent"};

  Matrix<T> B(*this);
  for (size_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] += A.vals_[i];
  return B;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) const {
  if (ncols_ != A.get_ncols() || nrows_ != A.get_nrows())
    throw std::runtime_error{"error: mult_mat: Matrices are not consistent"};

  Matrix<T> B(*this);
  for (size_t i = 0; i < vals_.size(); ++i)
    B.vals_[i] -= A.vals_[i];
  return B;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& A) {
  if (ncols_ != A.get_ncols() || nrows_ != A.get_nrows())
    throw std::runtime_error{"error: mult_mat: Matrices are not consistent"};

  for (size_t i = 0; i < vals_.size(); ++i)
    vals_[i] += A.vals_[i];
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& A) {
  if (ncols_ != A.get_ncols() || nrows_ != A.get_nrows())
    throw std::runtime_error{"error: mult_mat: Matrices are not consistent"};

  for (size_t i = 0; i < vals_.size(); ++i)
    vals_[i] -= A.vals_[i];
  return *this;
}

template <typename T>
size_t Matrix<T>::get_ncols() const { return ncols_; }

template <typename T>
size_t Matrix<T>::get_nrows() const { return nrows_; }

template <typename T>
std::vector<T> Matrix<T>::get_vals() const { return vals_; }

template <typename T1>
std::ostream& operator<<(std::ostream& out, const Matrix<T1>& A) {
  for (size_t i = 0; i < A.nrows_; ++i) {
    for (size_t j = 0; j < A.ncols_; ++j)
      out << A.vals_[i + j * A.nrows_] << "\t";
    out << "\n";
  }
  return out;
}

template <typename T>
std::istream& operator>>(std::istream& in, const Matrix<T>& A) {
  for (size_t j = 0; j < A.ncols_; ++j)
    for (size_t i = 0; i < A.nrows_; ++i)
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

template <typename T1>
Matrix<T1> operator*(const Matrix<T1>& A, const Matrix<T1>& B) {
  if (A.ncols_ != B.nrows_)
    throw std::runtime_error{"error: mult_mat: Matrices are not consistent"};

  Matrix<T1> C(A.nrows_, B.ncols_);

  // column-major order: A(i, j) = A.val[i + j * nrows_]
  for(size_t j = 0; j < B.ncols_; ++j)
    for(size_t k = 0; k < A.ncols_; ++k)
      for(size_t i = 0; i < A.nrows_; ++i)
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