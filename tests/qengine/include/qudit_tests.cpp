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

#include "qudit.h"
#include "types.h"

class QuditTests : public ::testing::Test {};

TEST_F(QuditTests, probability) {
  qengine::Qudit<double> qudit(2);

  EXPECT_EQ(qudit.get_probability(0), 1.0);
}

TEST_F(QuditTests, probabilities) {
  qengine::Qudit<double> qudit(3);

  EXPECT_EQ(qudit.get_probabilities(), qengine::RVec<double>({1.0, 0.0, 0.0}));
}

TEST_F(QuditTests, dagger) {
  qengine::Qudit<double> a(2);
  qengine::Qudit<double> b(2);

  EXPECT_EQ(a.dagger(), b);
}

TEST_F(QuditTests, braket) {
  qengine::Qudit<double> a(2);
  qengine::Qudit<double> b(2);

  EXPECT_EQ(b.times(a), 1.0);
}

TEST_F(QuditTests, ketbra) {
  qengine::Qudit<double> a(2);
  qengine::Qudit<double> b(2);
  qengine::CMat<double> M(2, {1.0, 0.0, 0.0, 0.0});

  EXPECT_EQ(a.tensor_times(b), M);
}