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

#ifndef NEW_FUNCTIONS_H
#define NEW_FUNCTIONS_H

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/value_vector.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/geometry/grid_uniform_quad_cache.h>

IGA_NAMESPACE_OPEN

template <int, int, int> class FunctionElement;

template<int dim, int range = 1, int rank = 1>
class NewFunction : public GridUniformQuadCache<dim>
{
public:
    using ElementAccessor = FunctionElement<dim, range, rank>;
    using ElementIterator = GridForwardIterator<ElementAccessor>;

    /** Types for the input/output evaluation arguments */
    ///@{
    /**
     * Type for the input argument of the function.
     */
    using Point = Points<dim>;

    /**
     * Type for the return of the function.
     */
    using Value = Values<dim, range, rank>;

    /**
     * Type for the derivative of the function.
     */
    template <int order>
    using Derivative = Derivatives<dim, range, rank, order>;

    /**
     * Type for the gradient of the function.
     */
    using Gradient = Derivative<1>;

    /**
     * Type for the hessian of the function.
     */
    using Hessian = Derivative<2>;

    /**
     * Type for the divergence of function.
     */
    using Div = Values<dim, range, rank-1>;
    ///@}

    /** @name Constructors and destructor. */
    ///@{
    /** Constructor */
    NewFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
                const ValueFlags flag,
                const Quadrature<dim> &quad);

    /** Destructor */
    virtual ~NewFunction() = default;
    ///@}


    virtual void init_element(ElementIterator &elem) = 0;
    virtual void fill_element(ElementIterator &elem) = 0;

protected:
    std::shared_ptr<typename ElementAccessor::CacheType>
    &get_cache(ElementIterator &elem);
};

IGA_NAMESPACE_CLOSE

#endif
