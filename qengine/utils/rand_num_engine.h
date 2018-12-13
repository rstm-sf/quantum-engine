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

#ifndef QENGINE_UTILS_RAND_NUM_ENGINE_H_
#define QENGINE_UTILS_RAND_NUM_ENGINE_H_

#include <random>

class RandNumEngine {
public:
  ~RandNumEngine();
  RandNumEngine(const RandNumEngine&) = delete;
  RandNumEngine(RandNumEngine&&) = delete;
  RandNumEngine& operator=(const RandNumEngine&) = delete;
  RandNumEngine& operator=(RandNumEngine&&) = delete;

  RandNumEngine(uint32_t seed = 1);

  std::mt19937& mte();

private:
  std::mt19937 mt_;
};

RandNumEngine::~RandNumEngine() = default;

RandNumEngine::RandNumEngine(uint32_t seed) : mt_{seed} {}

std::mt19937& RandNumEngine::mte() { return mt_; }

#endif // QENGINE_UTILS_RAND_NUM_ENGINE_H_