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