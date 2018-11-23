#include <complex>

#include <gtest/gtest.h>

#include "matrix.h"
#include "types.h"

class MatrixTests : public ::testing::Test {};

TEST_F(MatrixTests, simple) {
  qengine::VectorC vals({1.0, 0.0, 0.0, 1.0});
  qengine::MatrixC A(2, 2, vals);

  EXPECT_EQ(A.get_vals(), vals);
}

TEST_F(MatrixTests, scalar_product) {
  qengine::MatrixC A(2, 2, {1.0, 0.0, 0.0, 1.0});
  A = A * 10;

  EXPECT_EQ(A, qengine::MatrixC (2, 2, {10.0, 0.0, 0.0, 10.0}));
}

TEST_F(MatrixTests, minus) {
  qengine::MatrixC A(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::MatrixC B(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::MatrixC C(A - B);

  EXPECT_EQ(C, qengine::MatrixC (2, 2, {0.0, 0.0, 0.0, 0.0}));
}

TEST_F(MatrixTests, minus_equal) {
  qengine::MatrixC A(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::MatrixC B(2, 2, {1.0, 2.0, 3.0, 4.0});
  A -= B;

  EXPECT_EQ(A, qengine::MatrixC (2, 2, {0.0, 0.0, 0.0, 0.0}));
}

TEST_F(MatrixTests, plus) {
  qengine::MatrixC A(2, 2, {0.0, 0.0, 0.0, 0.0});
  qengine::MatrixC B(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::MatrixC C(A + B);

  EXPECT_EQ(C, qengine::MatrixC (2, 2, {1.0, 2.0, 3.0, 4.0}));
}

TEST_F(MatrixTests, plus_equal) {
  qengine::MatrixC A(2, 2, {0.0, 0.0, 0.0, 0.0});
  qengine::MatrixC B(2, 2, {1.0, 2.0, 3.0, 4.0});
  A += B;

  EXPECT_EQ(A, qengine::MatrixC (2, 2, {1.0, 2.0, 3.0, 4.0}));
}

TEST_F(MatrixTests, matrix_product) {
  qengine::MatrixC A(2, 2, {1.0, 2.0, 3.0, 4.0});
  qengine::MatrixC B(2, 2, {4.0, 3.0, 2.0, 1.0});
  qengine::MatrixC C(A * B);

  EXPECT_EQ(C, qengine::MatrixC (2, 2, {13.0, 20.0, 5.0, 8.0}));
}

TEST_F(MatrixTests, transpose) {
  qengine::MatrixC A(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});

  EXPECT_EQ(A.transpose(), qengine::MatrixC (3, 2, {1.0, 3.0, 5.0,
                                                    2.0, 4.0, 6.0 }));
}

TEST_F(MatrixTests, dagger) {
  qengine::MatrixC A(2, 3, {1.0 + 1.0i, 2.0 + 1.0i,
                            3.0 + 1.0i, 4.0 + 1.0i,
                            5.0 + 1.0i, 6.0 + 1.0i});

  EXPECT_EQ(A.dagger(), qengine::MatrixC (
    3, 2, {1.0 - 1.0i, 3.0 - 1.0i, 5.0 - 1.0i,
           2.0 - 1.0i, 4.0 - 1.0i, 6.0 - 1.0i }));
}

TEST_F(MatrixTests, conjugate) {
  qengine::MatrixC A(2, 3, {1.0 + 1.0i, 2.0 + 1.0i,
                            3.0 + 1.0i, 4.0 + 1.0i,
                            5.0 + 1.0i, 6.0 + 1.0i});

  EXPECT_EQ(A.conjugate(), qengine::MatrixC (2, 3, {1.0 - 1.0i, 2.0 - 1.0i,
                                                    3.0 - 1.0i, 4.0 - 1.0i,
                                                    5.0 - 1.0i, 6.0 - 1.0i}));
}

TEST_F(MatrixTests, trace) {
  qengine::MatrixC A(2, 2, {1.0 + 1.0i, 2.0 + 1.0i,
                            3.0 + 1.0i, 4.0 + 1.0i});

  EXPECT_EQ(A.trace(), static_cast<qengine::Complex>(5.0 + 2.0i));
}

TEST_F(MatrixTests, tensor_product) {
  qengine::MatrixC A(2, 2, {1.0, 2.0,
                            3.0, 4.0 });
  qengine::MatrixC B(2, 2, {4.0, 3.0,
                            2.0, 1.0 });
  qengine::MatrixC C(4, 4, { 4.0, 3.0,  8.0,  6.0,
                             2.0, 1.0,  4.0,  2.0,
                            12.0, 9.0, 16.0, 12.0,
                             6.0, 3.0,  8.0,  4.0 });

  EXPECT_EQ(A.tensor_times(B), C);
}