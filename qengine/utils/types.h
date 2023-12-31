// Copyright (C) 2018 - 2019 Rustam Sayfutdinov (rstm.sf@gmail.com)
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

#ifndef QENGINE_UTILS_TYPE_H_
#define QENGINE_UTILS_TYPE_H_

#include <complex>
#include <cstdint>
#include <vector>

#include "square_matrix.h"

namespace qengine {
inline namespace type {

using CReg = uint32_t;

template <typename T>
using RVec = std::vector<T>;
template <typename T>
using RMat = SquareMatrix<T>;

template <typename T>
using Cmplx = std::complex<T>;

template <typename T>
using CVec = std::vector<Cmplx<T>>;
template <typename T>
using CMat = SquareMatrix<Cmplx<T>>;

using FCmplx = Cmplx<float>;
using DCmplx = Cmplx<double>;

using FCVec = CVec<float>;
using DCVec = CVec<double>;

using FCMat = CMat<float>;
using DCMat = CMat<double>;

} // namespace type
} // namespace qengine

#endif // QENGINE_UTILS_TYPE_H_