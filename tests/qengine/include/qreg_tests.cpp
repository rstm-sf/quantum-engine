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

#include "qreg.h"
#include "types.h"

class QRegTests : public ::testing::Test {};

TEST_F(QRegTests, probabilities) {
  qengine::QReg<double> qudit(3);

  EXPECT_EQ(qudit.probabilities(), qengine::RVec<double>({1.0, 0.0, 0.0}));
}

TEST_F(QRegTests, conjugate) {
  qengine::QReg<double> a(2);
  qengine::QReg<double> b(2);

  EXPECT_EQ(a.conjugate().probabilities(), b.probabilities());
}

TEST_F(QRegTests, braket) {
  qengine::QReg<double> a(2);
  qengine::QReg<double> b(2);

  EXPECT_EQ(b.braket_product(a), 1.0);
}

TEST_F(QRegTests, ketbra) {
  qengine::QReg<double> a(2);
  qengine::QReg<double> b(2);
  qengine::CMat<double> M(2, {1.0, 0.0, 0.0, 0.0});

  EXPECT_EQ(a.ketbra_product(b), M);
}

TEST_F(QRegTests, apply) {
  qengine::QReg<double> a(3);
  qengine::CMat<double> M(3, { 0.0, 0.0, 1.0,
                               0.0, 1.0, 0.0,
                               1.0, 0.0, 0.0 });
  std::vector<double> probs({0.0, 0.0, 1.0});
  a.apply(M);

  EXPECT_EQ(a.probabilities(), probs);
}