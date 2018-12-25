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

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "circuit.h"
#include "qreg.h"

class CircuitTests : public ::testing::Test {};

TEST_F(CircuitTests, simple) {
  uint64_t nreg = 2;
  uint64_t dim = 3;
  qengine::Circuit<double> circuit(nreg, dim);

  EXPECT_EQ(
    circuit.qregs()[0].probabilities(),
    qengine::QReg<double>(dim).probabilities());
}

TEST_F(CircuitTests, measure) {
  uint64_t nreg = 2;
  uint64_t dim = 3;
  qengine::Circuit<double> circuit(nreg, dim);

  qengine::RMat<double> op(dim, { 0.0, 0.0, 1.0,
                                  0.0, 1.0, 0.0,
                                  1.0, 0.0, 0.0 });
  circuit.apply(0, op);
  circuit.measure(0, 0);

  EXPECT_EQ(circuit.cregs()[0], 2);
}

TEST_F(CircuitTests, applyX) {
  uint64_t nreg = 2;
  uint64_t dim = 3;
  qengine::Circuit<double> circuit(nreg, dim);

  circuit.applyX(0);
  circuit.measure(0, 0);

  EXPECT_EQ(circuit.cregs()[0], 2);
}

TEST_F(CircuitTests, applyZ) {
  uint64_t nreg = 2;
  uint64_t dim = 2;
  qengine::Circuit<double> circuit(nreg, dim);

  circuit.applyX(0);
  circuit.applyZ(0);

  EXPECT_EQ(
    circuit.qregs()[0].probabilities(),
    std::vector<double>({0.0, 1.0}));
}

TEST_F(CircuitTests, applyF_1_F) {
  uint64_t nreg = 2;
  uint64_t dim = 2;
  qengine::Circuit<double> circuit(nreg, dim);
  std::vector<double> probs({1.0, 0.0, 0.0, 0.0});
  bool has_fail = false;

  circuit.applyF(0);
  circuit.applyF(0);

  for (std::size_t i = 0; i < 2; ++i)
    if (circuit.qregs()[0].probabilities()[i] - probs[i] > 1.0e-5) {
      has_fail = true;
      break;
    }

  EXPECT_TRUE(has_fail == false);
}

TEST_F(CircuitTests, applyF_1_F_dim4) {
  uint64_t nreg = 2;
  uint64_t dim = 4;
  qengine::Circuit<double> circuit(nreg, dim);
  std::vector<double> probs({1.0, 0.0, 0.0, 0.0});
  bool has_fail = false;

  circuit.applyF(0);
  circuit.applyFconjugate(0);

  for (std::size_t i = 0; i < dim; ++i)
    if (circuit.qregs()[0].probabilities()[i] - probs[i] > 1.0e-5) {
      has_fail = true;
      break;
    }

  EXPECT_TRUE(has_fail == false);
}