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

/*
 *  Test of
 *  author: pauletti
 *  date: 2014-08-18
 *
 */

#include "../tests.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_element_handler.h>

#include <boost/variant.hpp>



//class function {
//public:
//    template<int k>
//    virtual void reset() const;
//};


template<int dim>
class Function : public GridElementHandler<dim>
{
    using variant_type = boost::variant<
            Quadrature<1>,
            Quadrature<2> > ;


    using parent_t = GridElementHandler<dim>;
public:
    Function(std::shared_ptr<const CartesianGrid<dim>> grid)
    :
        GridElementHandler<dim>(grid)
        {}


    virtual void reset(const NewValueFlags &flag, const variant_type& quad)
    {
        impl.flag = flag;
        impl.grid = this;
        boost::apply_visitor(impl, quad);
    }

private:
    // A class that works with boost::apply_visitor
    struct dispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            grid->template reset<T::dim>(flag, quad);
            quad.print_info(out);
        }

        NewValueFlags flag;
        parent_t *grid;
    };
    dispatcher impl;
};


int main() {
    const int dim = 2;
    auto grid = CartesianGrid<dim>::create(3);
    Function<2> x(grid);

    x.reset(NewValueFlags::none, QGauss<1>(2));
    x.reset(NewValueFlags::none, QGauss<2>(2));

}

