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

#include "math_operations.h"
#include "types.h"

namespace qengine {
inline namespace qstate {

template <typename T>
class Qudit {
public:
  Qudit<T>();
  ~Qudit<T>();
  Qudit<T>(const Qudit<T>&);
  Qudit<T>(Qudit<T>&&);
  Qudit<T>& operator=(const Qudit<T>&);
  Qudit<T>& operator=(Qudit<T>&&);

  Qudit<T>(uint64_t dim, State state = State::KET);

  uint64_t get_dim() const;
  CVec<T> get_amplitudes() const;
  State get_state() const;

  T get_probability(uint64_t i) const;
  RVec<T> get_probabilities() const;

  Qudit<T> dagger() const;
  Cmplx<T> times(const Qudit<T>& ket) const;
  CMat<T> tensor_times(const Qudit<T>& ket) const;

  template <typename S>
  void apply(const Matrix<S>& mat);

  template <typename S>
  friend bool operator==(const Qudit<S>& a, const Qudit<S>& b);

  template <typename S>
  friend bool operator!=(const Qudit<S>& a, const Qudit<S>& b);

private:
  uint64_t dim_;
  CVec<T> amplitudes_;
  State state_;
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
Qudit<T>::Qudit(uint64_t dim, State state)
  : dim_{dim}, state_{state}, amplitudes_(dim) {
  amplitudes_[0] = 1.0;
}

template <typename T>
uint64_t Qudit<T>::get_dim() const { return dim_; }

template <typename T>
CVec<T> Qudit<T>::get_amplitudes() const { return amplitudes_; }

template <typename T>
State Qudit<T>::get_state() const { return state_; }

template <typename T>
T Qudit<T>::get_probability(uint64_t i) const {
  return std::real(amplitudes_[i] * std::conj(amplitudes_[i]));
}

template <typename T>
RVec<T> Qudit<T>::get_probabilities() const {
  RVec<T> probabilities;
  probabilities.reserve(dim_);
  for (uint64_t i = 0; i < dim_; ++i)
    probabilities.push_back(get_probability(i));
  return probabilities;
}

template <typename T>
Qudit<T> Qudit<T>::dagger() const {
  Qudit<T> res(*this);
  res.state_ = (state_ == State::KET) ? State::BRA : State::KET;
  for (auto & res_a : res.amplitudes_)
    res_a = std::conj(res_a);
  return res;
}

template <typename T>
Cmplx<T> Qudit<T>::times(const Qudit<T>& ket) const {
  Expects(state_ == State::BRA && ket.state_ == State::KET);
  Expects(ket.dim_ == dim_);

  Cmplx<T> res(0.0);
  for(uint64_t i = 0; i < dim_; ++i)
    res += amplitudes_[i] * ket.amplitudes_[i];
  return res;
}

template <typename T>
CMat<T> Qudit<T>::tensor_times(const Qudit<T>& bra) const {
  Expects(state_ == State::KET && bra.state_ == State::BRA);

  return ketbra_tensor_product(amplitudes_, bra.amplitudes_);
}

template <typename T>
template <typename S>
void Qudit<T>::apply(const Matrix<S>& mat) {
  amplitudes_ = (state_ == State::KET) ? mat * amplitudes_ : amplitudes_ * mat;
}

template <typename S>
bool operator==(const Qudit<S>& a, const Qudit<S>& b) {
  return (a.state_ == b.state_) && (a.amplitudes_ == b.amplitudes_);
}

template <typename S>
bool operator!=(const Qudit<S>& a, const Qudit<S>& b) { return !(a == b); }

} // namespace qstate
} // namespace qengine

#endif // QENGINE_INCLUDE_QUDIT_H_