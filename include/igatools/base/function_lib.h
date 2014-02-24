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

#ifndef __IGA_FUNCTION_LIB_H_
#define __IGA_FUNCTION_LIB_H_

#include <igatools/base/config.h>
#include <igatools/base/function.h>

IGA_NAMESPACE_OPEN

/**
 * Collection of useful functions derived from the Function class.
 */
namespace functions
{

/**
 * Constant scalar function.
 */
template<int dim, int range=1, int rank=1>
class ConstantFunction
    : public Function<dim, range, rank>
{
public:
    using PointType = typename Function<dim, range, rank>::PointType;
    using ValueType = typename Function<dim, range, rank>::ValueType;

    /**
     * Construct a constant function with the given value.
     */
    ConstantFunction(const ValueType value);

    /**
     * Destructor.
     */
    virtual ~ConstantFunction();

    /**
     * Compute the @p values of Function at some @p points.
     */
    void
    evaluate(const std::vector<PointType> &points,
             std::vector<ValueType> &values) const;

private:
    /** Constant given value. */
    const ValueType value_;

};

} // of namespace functions.


IGA_NAMESPACE_CLOSE

#endif /* __IGA_FUNCTION_LIB_H_ */
