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

#include <gtest/gtest.h>

#include "math_operations.h"
#include "matrix.h"
#include "types.h"

class MathOperationsTests : public ::testing::Test {};

TEST_F(MathOperationsTests, ketbra_tensor_product) {
  qengine::VectorC ket {{1.0, 0.0, 1.0}};
  qengine::VectorC bra {{3.0, 2.0, 1.0}};
  qengine::MatrixC A = qengine::ketbra_tensor_product(ket, bra);
  qengine::MatrixC B(3, 3, {
    3.0, 0.0, 3.0,
    2.0, 0.0, 2.0,
    1.0, 0.0, 1.0
  });

  EXPECT_EQ(A, B);
}