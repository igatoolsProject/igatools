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

#ifndef IG_FUNCTIONS_H
#define IG_FUNCTIONS_H

#include <igatools/base/new_function.h>
#include <igatools/linear_algebra/distributed_vector.h>

IGA_NAMESPACE_OPEN
template<class Space>
class IgFunction : public NewFunction<Space::dim, Space::codim, Space::range, Space::rank>
{
public:
    static const int dim = Space::dim;
    static const int codim = Space::codim;
    static const int range = Space::range;
    static const int rank = Space::rank;

private:
    using parent_t = NewFunction<dim, codim, range, rank>;

public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    using CoeffType = Vector<LAPack::trilinos>;

    IgFunction(const NewValueFlags &flag, const Quadrature<dim> &quad,
               std::shared_ptr<const Space> space,
               const CoeffType &coeff);

    void init_elem(ElementAccessor &elem);

    void fill_elem(ElementAccessor &elem);

private:
    FunctionFlags flag_;

    Quadrature<dim> quad_;

    std::shared_ptr<const Space> space_;

    const CoeffType coeff_;

    typename Space::ElementIterator elem_;

    typename Space::ElementHandler space_filler_;
};

IGA_NAMESPACE_CLOSE

#endif
