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

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

#include "circuit.h"
#include "qreg.h"

int main(int argc, char const *argv[]) {
  uint64_t n = 8;
  float epsilon = 0.4f;
  uint64_t dim = 32;

  std::vector<uint64_t> K({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
  uint64_t word_a = 16;
  uint64_t word_b = 16;

  uint64_t nreg = 1;
  qengine::Circuit<double> circuit(nreg, dim);

  std::cout << "n   = " << n << "\n";
  std::cout << "dim = " << dim << "\n";
  std::cout << "a   = " << word_a << "\n";
  std::cout << "b   = " << word_b << "\n";
  std::cout << "K   = {";
  for (const auto& k : K)
    std::cout << k << ", ";
  std::cout << "}" << std::endl;

  // Применение хеш-функции
  // Получаем состяние $\ket{\psi(a)}$
  circuit.applyF(0);
  for (const auto& k : K)
    circuit.applyZ(0, word_a * k);


  // Reverse-тест
  std::reverse(std::begin(K), std::end(K));
  for (const auto& k : K)
    circuit.applyZ(0, word_b * k);
  circuit.applyFconjugate(0);
  circuit.measure(0, 0);

  auto creg = circuit.cregs();

  if (creg[0] == 0LL) {
    std::cout << "\nTEST PASSED" << std::endl;
  } else {
    std::cout << "\nTEST FAIL" << std::endl;
  }

  return 0;
}