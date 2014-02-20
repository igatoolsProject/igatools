//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
//
// This file is part of the igatools library.
//
// The igatools library is free software: you can use it, redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-+--------------------------------------------------------------------

#include <igatools/base/tensor.h>
#include <igatools/base/tensor-inline.h>

IGA_NAMESPACE_OPEN

Tdouble::
Tdouble(const Real val)
{
    Assert(!std::isnan(val),ExcNotANumber());
    Assert(!std::isinf(val),ExcNumberNotFinite());

    val_ = val;
}



Tdouble &
Tdouble::
operator=(const value_t &val)
{
    Assert(!std::isnan(val),ExcNotANumber());
    Assert(!std::isinf(val),ExcNumberNotFinite());

    val_ = val;
    return (*this);
}



auto
Tdouble::
operator[](const int i) noexcept -> value_t &
{
    return *this;
}



auto
Tdouble::
operator[](const int i) const noexcept -> const value_t &
{
    return *this;
}



auto
Tdouble::
operator()(const product_Index  &i) noexcept -> value_t &
{
    return *this;
}



auto
Tdouble::
operator()(const product_Index  &i) const noexcept -> const value_t &
{
    return *this;
}



auto
Tdouble::
operator()(const int i) noexcept -> value_t &
{
    return *this;
}



auto
Tdouble::
operator()(const int i) const noexcept -> const value_t &
{
    return *this;
}



Tdouble &
Tdouble::
operator+=(const Real val) noexcept
{
    Assert(!std::isnan(val),ExcNotANumber());
    Assert(!std::isinf(val),ExcNumberNotFinite());

    val_ += val;
    return (*this);
}



Tdouble &
Tdouble::
operator-=(const Real val) noexcept
{
    Assert(!std::isnan(val),ExcNotANumber());
    Assert(!std::isinf(val),ExcNumberNotFinite());

    val_ -= val;
    return (*this);
}



Tdouble &
Tdouble::
operator*=(const Real val) noexcept
{
    Assert(!std::isnan(val),ExcNotANumber());
    Assert(!std::isinf(val),ExcNumberNotFinite());

    val_ *= val;
    return (*this);
}



Tdouble &
Tdouble::
operator/=(const Real val) noexcept
{
    Assert(!std::isnan(val),ExcNotANumber());
    Assert(!std::isinf(val),ExcNumberNotFinite());
    Assert(val != Real(0.0),ExcDivideByZero());

    val_ /= val;
    return (*this);
}



Real
Tdouble::
norm() const noexcept
{
    return std::abs(val_);
}



Real
Tdouble::
norm_square() const noexcept
{
    return val_ * val_;
}



auto
Tdouble::
flat_to_tensor_index(const int flat_index) const noexcept -> product_Index
{
    return product_Index();
}



int
Tdouble::
tensor_to_flat_index(const product_Index &tensor_index) const noexcept
{
    return 0;
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/tensor.inst>
