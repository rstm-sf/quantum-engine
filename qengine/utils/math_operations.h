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