// Copyright (C) 2018 Rustam Sayfutdinov, rstm.sf@gmail.com
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

#include "types.h"

namespace qengine {
inline namespace mo {

MatrixC get_SWAP() {
  // column-major order: A(i, j) = A.val[i + j * nrows]
  return MatrixC(2, 2, { 1.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 1.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 1.0 });
}

MatrixC get_I(size_t n) {
  if (n == 1) {
    return MatrixC(1, 1, {1.0});
  } else if (n == 2) {
    return MatrixC(2, 2, {1.0, 0.0, 0.0, 1.0});
  } else {
    MatrixC I(n, n);
    for (size_t i = 0; i < n; ++i)
        I(i, i) = 1.0;
    return I;
  }
}

void identity(MatrixC & I) {
  Expects(I.get_nrows() == I.get_ncols());

  I = get_I(I.get_nrows());
}

MatrixC ketbra_tensor_product(const VectorC& ket, const VectorC& bra) {
  MatrixC res(ket.size(), bra.size());
  for (size_t j = 0; j < bra.size(); ++j)
    for (size_t i = 0; i < ket.size(); ++i)
      res(i, j) = ket[i] * bra[j];
  return res;
}

} // namespace mo
} // namespace qengine

#endif // QENGINE_UTILS_MATH_OPERATIONS_H_