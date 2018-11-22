#ifndef QENGINE_UTILS_TYPE_H_
#define QENGINE_UTILS_TYPE_H_

#include <complex>
#include <vector>

#include "matrix.h"

namespace qengine {
inline namespace type {

using Complex = std::complex<double>;
using VectorC = std::vector<Complex>;
using MatrixC = Matrix<Complex>;

} // namespace type
} // namespace qengine

#endif // QENGINE_UTILS_TYPE_H_