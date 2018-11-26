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
  using DCVec = std::vector<std::complex<double>>;
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCVec vals({1.0, 0.0, 0.0, 1.0});
  DCMat A(2, 2, vals);

  EXPECT_EQ(A.get_vals(), vals);
}

TEST_F(MatrixTests, scalar_product) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 2, {1.0, 0.0, 0.0, 1.0});
  A = A * 10;

  EXPECT_EQ(A, DCMat (2, 2, {10.0, 0.0, 0.0, 10.0}));
}

TEST_F(MatrixTests, minus) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 2, {1.0, 2.0, 3.0, 4.0});
  DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  DCMat C(A - B);

  EXPECT_EQ(C, DCMat (2, 2, {0.0, 0.0, 0.0, 0.0}));
}

TEST_F(MatrixTests, minus_equal) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 2, {1.0, 2.0, 3.0, 4.0});
  DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  A -= B;

  EXPECT_EQ(A, DCMat (2, 2, {0.0, 0.0, 0.0, 0.0}));
}

TEST_F(MatrixTests, plus) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 2, {0.0, 0.0, 0.0, 0.0});
  DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  DCMat C(A + B);

  EXPECT_EQ(C, DCMat (2, 2, {1.0, 2.0, 3.0, 4.0}));
}

TEST_F(MatrixTests, plus_equal) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 2, {0.0, 0.0, 0.0, 0.0});
  DCMat B(2, 2, {1.0, 2.0, 3.0, 4.0});
  A += B;

  EXPECT_EQ(A, DCMat (2, 2, {1.0, 2.0, 3.0, 4.0}));
}

TEST_F(MatrixTests, matrix_product) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(3, 5, { 1.0,  6.0, 11.0,
                  2.0,  7.0, 12.0,
                  3.0,  8.0, 13.0,
                  4.0,  9.0, 14.0,
                  5.0, 10.0, 15.0 });
  DCMat B(5, 3, { 15.0, 12.0, 9.0, 6.0, 3.0,
                  14.0, 11.0, 8.0, 5.0, 2.0,
                  13.0, 10.0, 7.0, 4.0, 1.0 });
  DCMat C(3, 3, { 105.0, 330.0, 555.0,
                   90.0, 290.0, 490.0,
                   75.0, 250.0, 425.0 });

  EXPECT_EQ(A * B, C);
}

TEST_F(MatrixTests, transpose) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});

  EXPECT_EQ(A.transpose(), DCMat (3, 2, { 1.0, 3.0, 5.0,
                                          2.0, 4.0, 6.0 }));
}

TEST_F(MatrixTests, dagger) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 3, { 1.0 + 1.0i, 2.0 + 1.0i,
                  3.0 + 1.0i, 4.0 + 1.0i,
                  5.0 + 1.0i, 6.0 + 1.0i});

  EXPECT_EQ(A.dagger(), DCMat (
    3, 2, {1.0 - 1.0i, 3.0 - 1.0i, 5.0 - 1.0i,
           2.0 - 1.0i, 4.0 - 1.0i, 6.0 - 1.0i }));
}

TEST_F(MatrixTests, conjugate) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 3, { 1.0 + 1.0i, 2.0 + 1.0i,
                  3.0 + 1.0i, 4.0 + 1.0i,
                  5.0 + 1.0i, 6.0 + 1.0i});

  EXPECT_EQ(A.conjugate(), DCMat (2, 3, {1.0 - 1.0i, 2.0 - 1.0i,
                                         3.0 - 1.0i, 4.0 - 1.0i,
                                         5.0 - 1.0i, 6.0 - 1.0i}));
}

TEST_F(MatrixTests, trace) {
  using DCmplx = std::complex<double>;
  using DCMat = qengine::Matrix<DCmplx>;

  DCMat A(2, 2, {1.0 + 1.0i, 2.0 + 1.0i,
                 3.0 + 1.0i, 4.0 + 1.0i});

  EXPECT_EQ(A.trace(), static_cast<DCmplx>(5.0 + 2.0i));
}

TEST_F(MatrixTests, tensor_product) {
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(2, 2, { 1.0, 2.0,
                  3.0, 4.0 });
  DCMat B(2, 2, { 4.0, 3.0,
                  2.0, 1.0 });
  DCMat C(4, 4, { 4.0, 3.0,  8.0,  6.0,
                  2.0, 1.0,  4.0,  2.0,
                 12.0, 9.0, 16.0, 12.0,
                  6.0, 3.0,  8.0,  4.0 });

  EXPECT_EQ(A.tensor_times(B), C);
}

TEST_F(MatrixTests, mat_vec_product_r) {
  using DCVec = std::vector<std::complex<double>>;
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(4, 4, { 1.0, 5.0,  9.0, 13.0,
                  2.0, 6.0, 10.0, 14.0,
                  3.0, 7.0, 11.0, 15.0,
                  4.0, 8.0, 12.0, 16.0 });
  DCVec b({4.0, 3.0, 2.0, 1.0});
  DCVec c({20.0, 60.0, 100.0, 140.0});

  EXPECT_EQ(A * b, c);
}

TEST_F(MatrixTests, mat_vec_product_l) {
  using DCVec = std::vector<std::complex<double>>;
  using DCMat = qengine::Matrix<std::complex<double>>;

  DCMat A(4, 4, { 1.0, 5.0,  9.0, 13.0,
                  2.0, 6.0, 10.0, 14.0,
                  3.0, 7.0, 11.0, 15.0,
                  4.0, 8.0, 12.0, 16.0 });
  DCVec b({4.0, 3.0, 2.0, 1.0});
  DCVec c({50.0, 60.0, 70.0, 80.0});

  EXPECT_EQ(b * A, c);
}