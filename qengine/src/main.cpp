#include <iostream>

#include "matrix.h"
#include "types.h"

int main(int /*argc*/, char* /*argv*/[]) {
  std::cout << "Hello!" << std::endl;
  qengine::MatrixC A(2, 2, {1.0, 2.0,
                            3.0, 4.0 });
  qengine::MatrixC B(2, 2, {4.0, 3.0,
                            2.0, 1.0 });
  std::cout << A.tensor_times(B) << std::endl;
  return 0;
}