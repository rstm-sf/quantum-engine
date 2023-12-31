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

  uint64_t nreg() const;
  uint64_t dim() const;

  void apply(uint64_t idx_qreg, RMat<T>);
  void apply(uint64_t idx_qreg, CMat<T>);

  void applyX(uint64_t idx_qreg, uint64_t i, T x, T y);
  void applyX(uint64_t idx_qreg, uint64_t i, Cmplx<T> x, Cmplx<T> y);
  void applyZ(uint64_t idx_qreg, uint64_t i, double tau);
  void applyXconjugate(uint64_t idx_qreg, uint64_t i, T x, T y);
  void applyXconjugate(uint64_t idx_qreg, uint64_t i, Cmplx<T> x, Cmplx<T> y);
  void applyZconjugate(uint64_t idx_qreg, uint64_t i, double tau);

  void measure(uint64_t idx_qreg, uint64_t idx_creg);

  CReg get_creg(uint64_t idx_creg) const;

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
uint64_t Circuit<T>::nreg() const { return qregs_.size(); }

template <typename T>
uint64_t Circuit<T>::dim() const { return qregs_[0].sdim(); }

template <typename T>
void Circuit<T>::apply(uint64_t idx_qreg, RMat<T> mat_op) {
  qregs_[idx_qreg].apply(mat_op);
}

template <typename T>
void Circuit<T>::apply(uint64_t idx_qreg, CMat<T> mat_op) {
  qregs_[idx_qreg].apply(mat_op);
}

template <typename T>
void Circuit<T>::applyX(uint64_t idx_qreg, uint64_t i, T x, T y) {
  qregs_[idx_qreg].applyX(i, x, y);
}

template <typename T>
void Circuit<T>::applyX(uint64_t idx_qreg, uint64_t i, Cmplx<T> x, Cmplx<T> y) {
  qregs_[idx_qreg].applyX(i, x, y);
}

template <typename T>
void Circuit<T>::applyZ(uint64_t idx_qreg, uint64_t i, double tau) {
  qregs_[idx_qreg].applyZ(i, tau);
}

template <typename T>
void Circuit<T>::applyXconjugate(
    uint64_t idx_qreg, uint64_t i, T x, T y) {
  qregs_[idx_qreg].applyXconjugate(i, x, y);
}

template <typename T>
void Circuit<T>::applyXconjugate(
    uint64_t idx_qreg, uint64_t i, Cmplx<T> x, Cmplx<T> y) {
  qregs_[idx_qreg].applyXconjugate(i, x, y);
}

template <typename T>
void Circuit<T>::applyZconjugate(uint64_t idx_qreg, uint64_t i, double tau) {
  qregs_[idx_qreg].applyZconjugate(i, tau);
}

template <typename T>
void Circuit<T>::measure(uint64_t idx_qreg, uint64_t idx_creg) {
  RVec<T> probs(qregs_[idx_qreg].probabilities());
  std::discrete_distribution<CReg> distrib(probs.begin(), probs.end());

  const CReg result = distrib(rand_eng_.mte());

  // TODO: применить ли изменение состояния после измерения?
  // RMat<T> Z(probs.size());
  // Z(result, result) = 1.0 / std::sqrt(probs[result]);
  // qregs_[idx_qreg].apply(Z);

  cregs_[idx_creg] = result;
}

template <typename T>
CReg Circuit<T>::get_creg(uint64_t idx_creg) const { return cregs_[idx_creg]; }

} // namespace qsystem
} // namespace qengine

#endif // QENGINE_INCLUDE_CIRCUIT_H_
