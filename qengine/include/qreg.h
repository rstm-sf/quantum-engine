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

#ifndef QENGINE_INCLUDE_QREG_H_
#define QENGINE_INCLUDE_QREG_H_

#include "ireg.h"
#include "math_operations.h"
#include "types.h"

namespace qengine {
inline namespace qstate {

template <typename T>
class QReg : public IReg<T> {
public:
  QReg<T>();
  virtual ~QReg<T>();
  QReg<T>(const QReg<T>&);
  QReg<T>(QReg<T>&&);
  QReg<T>& operator=(const QReg<T>&);
  QReg<T>& operator=(QReg<T>&&);

  QReg<T>(uint64_t sdim, uint64_t size = 1);

  virtual uint64_t size() const override;

  virtual void apply(const CMat<T> mat, uint64_t idx_qudit = 0);

  uint64_t dim() const;
  uint64_t sdim() const;

  RVec<T> probabilities() const;

  QReg<T> conjugate() const;
  Cmplx<T> braket_product(const QReg<T>& ket) const;
  CMat<T> ketbra_product(const QReg<T>& ket) const;

private:
  uint64_t sdim_;
  uint64_t size_;
  CVec<T> amplitudes_;
};

template <typename T>
QReg<T>::QReg() = default;

template <typename T>
QReg<T>::~QReg() = default;

template <typename T>
QReg<T>::QReg(const QReg<T>&) = default;

template <typename T>
QReg<T>::QReg(QReg<T>&&) = default;

template <typename T>
QReg<T>& QReg<T>::operator=(const QReg<T>&) = default;

template <typename T>
QReg<T>& QReg<T>::operator=(QReg<T>&&) = default;

template <typename T>
QReg<T>::QReg(uint64_t sdim, uint64_t size)
  : sdim_{sdim}, size_{size}, amplitudes_(sdim * size) {
  amplitudes_[0] = 1.0;
}

template <typename T>
uint64_t QReg<T>::size() const { return size_; }

template <typename T>
RVec<T> QReg<T>::probabilities() const {
  RVec<T> res;
  res.reserve(amplitudes_.size());
  for(auto const & a : amplitudes_)
    res.push_back(probability(a));
  return res;
}

template <typename T>
uint64_t QReg<T>::dim() const { return amplitudes_.size(); }

template <typename T>
uint64_t QReg<T>::sdim() const { return sdim_; }

template <typename T>
void QReg<T>::apply(const CMat<T> mat, uint64_t idx_qudit) {
  amplitudes_ = mat * amplitudes_;
}

template <typename T>
QReg<T> QReg<T>::conjugate() const {
  QReg<T> res(*this);
  for (auto & res_a : res.amplitudes_)
    res_a = std::conj(res_a);
  return res;
}

template <typename T>
Cmplx<T> QReg<T>::braket_product(const QReg<T>& ket) const {
  Expects(ket.amplitudes_.size() == amplitudes_.size());

  Cmplx<T> res(0.0);
  for(uint64_t i = 0; i < amplitudes_.size(); ++i)
    res += amplitudes_[i] * ket.amplitudes_[i];
  return res;
}

template <typename T>
CMat<T> QReg<T>::ketbra_product(const QReg<T>& bra) const {
  return ketbra_tensor_product(amplitudes_, bra.amplitudes_);
}

} // namespace qstate
} // namespace qengine

#endif // QENGINE_INCLUDE_QREG_H_