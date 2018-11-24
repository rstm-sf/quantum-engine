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

#ifndef QENGINE_UTILS_MATH_OPERATIONS_H_
#define QENGINE_UTILS_MATH_OPERATIONS_H_

#include <vector>

#include "matrix.h"

namespace qengine {
inline namespace mo {

template <typename T>
constexpr Matrix<T> SWAP_mat() {
  // column-major order: A(i, j) = A.val[i + j * nrows]
  return Matrix<T>(4, 4, { 1.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 0.0,
                           0.0, 1.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.0 });
}

template <typename T>
constexpr Matrix<T> I_mat_1x1() { return Matrix<T>(1, 1, {1.0}); }
template <typename T>
constexpr Matrix<T> I_mat_2x2() {
  return Matrix<T>(2, 2, {1.0, 0.0, 0.0, 1.0});
}

template <typename T>
Matrix<T> I_mat(size_t n) {
  if (n == 1) {
    return I_mat_1x1<T>();
  } else if (n == 2) {
    return I_mat_2x2<T>();
  } else {
    Matrix<T> res(n, n);
    for (size_t i = 0; i < n; ++i)
      res(i, i) = 1.0;
    return res;
  }
}

template <typename T>
Matrix<T> ketbra_tensor_product(
    const std::vector<T>& ket, const std::vector<T>& bra) {
  Matrix<T> res(ket.size(), bra.size());
  for (size_t j = 0; j < bra.size(); ++j)
    for (size_t i = 0; i < ket.size(); ++i)
      res(i, j) = ket[i] * bra[j];
  return res;
}

} // namespace mo
} // namespace qengine

#endif // QENGINE_UTILS_MATH_OPERATIONS_H_