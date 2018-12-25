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

#ifndef QENGINE_INCLUDE_CIRCUIT_H_
#define QENGINE_INCLUDE_CIRCUIT_H_

#include <cmath>
#include <cstdint>
#include <vector>

#include "rand_num_engine.h"
#include "qreg.h"
#include "types.h"

namespace qengine {
inline namespace qsystem {

template <typename T>
class Circuit {
public:
  Circuit<T>() = delete;
  ~Circuit<T>();
  Circuit<T>(const Circuit<T>&) = delete;
  Circuit<T>(Circuit<T>&&) = delete;
  Circuit<T>& operator=(const Circuit<T>&) = delete;
  Circuit<T>& operator=(Circuit<T>&&) = delete;

  Circuit<T>(uint64_t nreg, uint64_t dim);

  std::vector<QReg<T>> qregs() const;
  std::vector<CReg> cregs() const;

  void apply(uint64_t idx_qreg, RMat<T>);
  void apply(uint64_t idx_qreg, CMat<T>);

  void applyX(uint64_t idx_qreg, uint64_t x = 1);
  void applyZ(uint64_t idx_qreg, uint64_t z = 1);
  void applyF(uint64_t idx_qreg);

  void measure(uint64_t idx_qreg, uint64_t idx_creg);

protected:
  std::vector<QReg<T>> qregs_;
  std::vector<CReg> cregs_;
  RandNumEngine rand_eng_;
};

template <typename T>
Circuit<T>::~Circuit() = default;

template <typename T>
Circuit<T>::Circuit(uint64_t nreg, uint64_t dim)
  : qregs_(nreg, QReg<T>(dim)), cregs_(nreg), rand_eng_() {}

template <typename T>
std::vector<QReg<T>> Circuit<T>::qregs() const { return qregs_; }

template <typename T>
std::vector<CReg> Circuit<T>::cregs() const { return cregs_; }

template <typename T>
void Circuit<T>::apply(uint64_t idx_qreg, RMat<T> mat_op) {
  qregs_[idx_qreg].apply(mat_op);
}

template <typename T>
void Circuit<T>::apply(uint64_t idx_qreg, CMat<T> mat_op) {
  qregs_[idx_qreg].apply(mat_op);
}

template <typename T>
void Circuit<T>::applyX(uint64_t idx_qreg, uint64_t x) {
  qregs_[idx_qreg].applyX(x);
}

template <typename T>
void Circuit<T>::applyZ(uint64_t idx_qreg, uint64_t z) {
  qregs_[idx_qreg].applyZ(z);
}

template <typename T>
void Circuit<T>::applyF(uint64_t idx_qreg) {
  qregs_[idx_qreg].applyF();
}

template <typename T>
void Circuit<T>::measure(uint64_t idx_qreg, uint64_t idx_creg) {
  RVec<T> probs(qregs_[idx_qreg].probabilities());
  std::discrete_distribution<uint16_t> distrib(probs.begin(), probs.end());

  const uint16_t result = distrib(rand_eng_.mte());

  // TODO: применить ли изменение состояния после измерения?
  // RMat<T> Z(probs.size());
  // Z(result, result) = 1.0 / std::sqrt(probs[result]);
  // qregs_[idx_qreg].apply(Z);

  cregs_[idx_creg] = result;
}

} // namespace qsystem
} // namespace qengine

#endif // QENGINE_INCLUDE_CIRCUIT_H_