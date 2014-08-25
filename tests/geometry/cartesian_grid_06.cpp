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
 *  Test for find algorith on grid container iterator
 *
 *  author: pauletti
 *  date:
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>


template<int dim>
bool IsOdd (const typename CartesianGrid<dim>::ElementAccessor &elem)
{
  return ((elem.get_flat_index()%2)==1);
}

template<int dim>
void do_test()
{
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
        n_knots(i) = 2*i+2;
    auto grid = CartesianGrid<dim>::create(n_knots);

    auto it = std::find_if (grid->begin(), grid->end(), IsOdd<dim>);

}


#include <type_traits>

template<class T>
using mytype = decltype(declval<T>().print_info(declval<LogStream &>()));

template<class T>
EnableIf<std::is_void<mytype<T>>::value, bool >
test(int)
{
	return true;
}

template<class T>
bool
test(long)
{
	return false;
}

class A
{
public:
	void print_info(LogStream & )
	{

	}
};


class B
{};

int main()
{
	out << test<A>(0) << endl;
	out << test<B>(0) << endl;

    return 0;
}
