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

#ifndef QENGINE_INCLUDE_IREG_H_
#define QENGINE_INCLUDE_IREG_H_

#include "types.h"

namespace qengine {
inline namespace qstate {

class IReg<T> {
public:
  virtual ~IReg<T>() = default;

  virtual uint64_t size() const = 0;

  virtual RVec<T> probabilities() const = 0;
  virtual T probability(uint64_t idx) const = 0;

  virtual IReg<T> conjugate() const = 0;
  virtual Cmplx<T> braket_product(const IReg<T>& ket) const = 0;
  virtual CMat<T> ketbra_product(const IReg<T>& bra) const = 0;
};

} // namespace qstate
} // namespace qengine

#endif // QENGINE_INCLUDE_IREG_H_