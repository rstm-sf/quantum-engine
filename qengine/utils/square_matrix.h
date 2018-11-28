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

#ifndef QENGINE_UTILS_SQUARE_MATRIX_H_
#define QENGINE_UTILS_SQUARE_MATRIX_H_

#include <cstdint>
#include <initializer_list>
#include <vector>

#include "matrix.h"

namespace qengine {
inline namespace util {

template <typename T>
class SquareMatrix : public Matrix<T> {
public:
  SquareMatrix<T>();
  virtual ~SquareMatrix<T>();
  SquareMatrix<T>(const SquareMatrix<T>&);
  SquareMatrix<T>(SquareMatrix<T>&&);
  SquareMatrix<T>& operator=(const SquareMatrix<T>&);
  SquareMatrix<T>& operator=(SquareMatrix<T>&&);

  SquareMatrix<T>(uint64_t n);
  SquareMatrix<T>(uint64_t n, const std::vector<T>& vals);
  SquareMatrix<T>(uint64_t n, const std::initializer_list<T>& vals);
};

template <typename T>
SquareMatrix<T>::SquareMatrix() = default;

template <typename T>
SquareMatrix<T>::~SquareMatrix() = default;

template <typename T>
SquareMatrix<T>::SquareMatrix(const SquareMatrix<T>&) = default;

template <typename T>
SquareMatrix<T>::SquareMatrix(SquareMatrix<T>&&) = default;

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator=(const SquareMatrix<T>&) = default;

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator=(SquareMatrix<T>&&) = default;

template <typename T>
SquareMatrix<T>::SquareMatrix(uint64_t n)
  : Matrix<T>(n, n) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(uint64_t n, const std::vector<T>& vals)
  : Matrix<T>(n, n, vals) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(uint64_t n, const std::initializer_list<T>& vals)
  : Matrix<T>(n, n, vals) {}

} // namespace util
} // namespace qengine

#endif // QENGINE_UTILS_SQUARE_MATRIX_H_