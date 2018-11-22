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