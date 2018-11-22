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