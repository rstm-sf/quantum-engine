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

#ifndef QENGINE_INCLUDE_CREG_H_
#define QENGINE_INCLUDE_CREG_H_

#include "ireg.h"
#include "types.h"

namespace qengine {
inline namespace cstate {

template <typename T = uint16_t>
class CReg : public IReg<T> {
public:
  CReg<T>();
  virtual ~CReg<T>();
  CReg<T>(const CReg<T>&);
  CReg<T>(CReg<T>&&);
  CReg<T>& operator=(const CReg<T>&);
  CReg<T>& operator=(CReg<T>&&);

  CReg<T>(uint64_t size);

  virtual uint64_t size() const override;

  RVec<T> vals() const;

  T& operator[](uint64_t idx);
  T operator[](uint64_t idx) const;

private:
  RVec<T> vals_;
};

template <typename T>
CReg<T>::CReg() = default;

template <typename T>
CReg<T>::~CReg() = default;

template <typename T>
CReg<T>::CReg(const CReg<T>&) = default;

template <typename T>
CReg<T>::CReg(CReg<T>&&) = default;

template <typename T>
CReg<T>& CReg<T>::operator=(const CReg<T>&) = default;

template <typename T>
CReg<T>& CReg<T>::operator=(CReg<T>&&) = default;

template <typename T>
CReg<T>::CReg(uint64_t size) : vals_(size) {}

template <typename T>
uint64_t CReg<T>::size() const { return vals_.size(); }

template <typename T>
RVec<T> CReg<T>::vals() const { return vals_; }

template <typename T>
T& CReg<T>::operator[](uint64_t idx) { return vals_[idx]; }

template <typename T>
T CReg<T>::operator[](uint64_t idx) const { return vals_[idx]; }

} // namespace cstate
} // namespace qengine

#endif // QENGINE_INCLUDE_CREG_H_