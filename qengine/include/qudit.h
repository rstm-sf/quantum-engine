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

#ifndef QENGINE_INCLUDE_QUDIT_H_
#define QENGINE_INCLUDE_QUDIT_H_

#include "ireg.h"
#include "math_operations.h"
#include "types.h"

namespace qengine {
inline namespace qstate {

template <typename T>
class Qudit : public IReg<T> {
public:
  Qudit<T>();
  ~Qudit<T>();
  Qudit<T>(const Qudit<T>&);
  Qudit<T>(Qudit<T>&&);
  Qudit<T>& operator=(const Qudit<T>&);
  Qudit<T>& operator=(Qudit<T>&&);

  Qudit<T>(uint64_t dim);

  virtual uint64_t dim() const;
  virtual uint64_t sdim() const;
  virtual uint64_t nstates() const;

  virtual RVec<T> probabilities() const;

  Qudit<T> conjugate() const;
  Cmplx<T> braket_product(const Qudit<T>& ket) const;
  CMat<T> ketbra_product(const Qudit<T>& ket) const;

  template <typename S>
  friend bool operator==(const Qudit<S>& a, const Qudit<S>& b);

  template <typename S>
  friend bool operator!=(const Qudit<S>& a, const Qudit<S>& b);

private:
  CVec<T> amplitudes_;
};

template <typename T>
Qudit<T>::Qudit() = default;

template <typename T>
Qudit<T>::~Qudit() = default;

template <typename T>
Qudit<T>::Qudit(const Qudit<T>&) = default;

template <typename T>
Qudit<T>::Qudit(Qudit<T>&&) = default;

template <typename T>
Qudit<T>& Qudit<T>::operator=(const Qudit<T>&) = default;

template <typename T>
Qudit<T>& Qudit<T>::operator=(Qudit<T>&&) = default;

template <typename T>
Qudit<T>::Qudit(uint64_t dim) : amplitudes_(dim) {
  amplitudes_[0] = 1.0;
}

template <typename T>
uint64_t Qudit<T>::dim() const { return amplitudes_.size(); }

template <typename T>
uint64_t Qudit<T>::sdim() const { return amplitudes_.size(); }

template <typename T>
uint64_t Qudit<T>::nstates() const { return 1; }

template <typename T>
RVec<T> Qudit<T>::probabilities() const {
  RVec<T> res;
  res.reserve(amplitudes_.size());
  for(auto const & a : amplitudes_)
    res.push_back(probability(a));
  return res;
}

template <typename T>
Qudit<T> Qudit<T>::conjugate() const {
  Qudit<T> res(*this);
  for (auto & res_a : res.amplitudes_)
    res_a = std::conj(res_a);
  return res;
}

template <typename T>
Cmplx<T> Qudit<T>::braket_product(const Qudit<T>& ket) const {
  Expects(ket.sdim() == this->sdim());

  Cmplx<T> res(0.0);
  for(uint64_t i = 0; i < amplitudes_.size(); ++i)
    res += amplitudes_[i] * ket.amplitudes_[i];
  return res;
}

template <typename T>
CMat<T> Qudit<T>::ketbra_product(const Qudit<T>& bra) const {
  return ketbra_tensor_product(amplitudes_, bra.amplitudes_);
}

template <typename S>
bool operator==(const Qudit<S>& a, const Qudit<S>& b) {
  return a.amplitudes_ == b.amplitudes_;
}

template <typename S>
bool operator!=(const Qudit<S>& a, const Qudit<S>& b) { return !(a == b); }

} // namespace qstate
} // namespace qengine

#endif // QENGINE_INCLUDE_QUDIT_H_