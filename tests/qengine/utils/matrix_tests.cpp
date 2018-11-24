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

#include <gtest/gtest.h>

#include "types.h"

class MatrixTests : public ::testing::Test {};

TEST_F(MatrixTests, simple) {
  qengine::DCVec vals({1.0, 0.0, 0.0, 1.0});
  qengine::DCMat A(2, 2, vals);

  EXPECT_EQ(A.get_vals(), vals);
}

TEST_F(MatrixTests, scalar_product) {
  qengine::DCMat A(2, 2, {1.0, 0.0, 0.0, 1.0});
  A = A * 10;

  EXPECT_EQ(A, qengine::DCMat (2, 2, {10.0, 0.0, 0.0, 10.0}));
}

TEST_F(MatrixTests, minus) {
  qengine::DCMat A(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::DCMat C(A - B);

  EXPECT_EQ(C, qengine::DCMat (2, 2, {0.0, 0.0, 0.0, 0.0}));
}

TEST_F(MatrixTests, minus_equal) {
  qengine::DCMat A(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  A -= B;

  EXPECT_EQ(A, qengine::DCMat (2, 2, {0.0, 0.0, 0.0, 0.0}));
}

TEST_F(MatrixTests, plus) {
  qengine::DCMat A(2, 2, {0.0, 0.0, 0.0, 0.0});
  qengine::DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::DCMat C(A + B);

  EXPECT_EQ(C, qengine::DCMat (2, 2, {1.0, 2.0, 3.0, 4.0}));
}

TEST_F(MatrixTests, plus_equal) {
  qengine::DCMat A(2, 2, {0.0, 0.0, 0.0, 0.0});
  qengine::DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  A += B;

  EXPECT_EQ(A, qengine::DCMat (2, 2, {1.0, 2.0, 3.0, 4.0}));
}

TEST_F(MatrixTests, matrix_product) {
  qengine::DCMat A(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::DCMat B(2, 2, {4.0, 3.0, 2.0, 1.0});
  qengine::DCMat C(A * B);

  EXPECT_EQ(C, qengine::DCMat (2, 2, {13.0, 20.0, 5.0, 8.0}));
}

TEST_F(MatrixTests, transpose) {
  qengine::DCMat A(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});

  EXPECT_EQ(A.transpose(), qengine::DCMat (3, 2, {1.0, 3.0, 5.0,
                                                  2.0, 4.0, 6.0 }));
}

TEST_F(MatrixTests, dagger) {
  qengine::DCMat A(2, 3, {1.0 + 1.0i, 2.0 + 1.0i,
                          3.0 + 1.0i, 4.0 + 1.0i,
                          5.0 + 1.0i, 6.0 + 1.0i});

  EXPECT_EQ(A.dagger(), qengine::DCMat (
    3, 2, {1.0 - 1.0i, 3.0 - 1.0i, 5.0 - 1.0i,
           2.0 - 1.0i, 4.0 - 1.0i, 6.0 - 1.0i }));
}

TEST_F(MatrixTests, conjugate) {
  qengine::DCMat A(2, 3, {1.0 + 1.0i, 2.0 + 1.0i,
                          3.0 + 1.0i, 4.0 + 1.0i,
                          5.0 + 1.0i, 6.0 + 1.0i});

  EXPECT_EQ(A.conjugate(), qengine::DCMat (2, 3, {1.0 - 1.0i, 2.0 - 1.0i,
                                                  3.0 - 1.0i, 4.0 - 1.0i,
                                                  5.0 - 1.0i, 6.0 - 1.0i}));
}

TEST_F(MatrixTests, trace) {
  qengine::DCMat A(2, 2, {1.0 + 1.0i, 2.0 + 1.0i,
                          3.0 + 1.0i, 4.0 + 1.0i});

  EXPECT_EQ(A.trace(), static_cast<qengine::DCmplx>(5.0 + 2.0i));
}

TEST_F(MatrixTests, tensor_product) {
  qengine::DCMat A(2, 2, {1.0, 2.0,
                          3.0, 4.0 });
  qengine::DCMat B(2, 2, {4.0, 3.0,
                          2.0, 1.0 });
  qengine::DCMat C(4, 4, { 4.0, 3.0,  8.0,  6.0,
                           2.0, 1.0,  4.0,  2.0,
                          12.0, 9.0, 16.0, 12.0,
                           6.0, 3.0,  8.0,  4.0 });

  EXPECT_EQ(A.tensor_times(B), C);
}