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

#include <igatools/base/function_lib.h>
#include <igatools/base/exceptions.h>

IGA_NAMESPACE_OPEN

using std::vector;

namespace functions
{
template<int dim, int range, int rank>
ConstantFunction<dim, range, rank>::
ConstantFunction(const Value value)
    :value_ {value}
{}



template<int dim, int range, int rank>
ConstantFunction<dim, range, rank>::
~ConstantFunction()
{}



template<int dim, int range, int rank>
void
ConstantFunction<dim, range, rank>::
evaluate(
    const vector<Point> &points,
    vector<Value> &values) const
{
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size())) ;

    for (auto &value : values)
        value = value_;
}

} // of namespace functions.


IGA_NAMESPACE_CLOSE

#include <igatools/base/function_lib.inst>
