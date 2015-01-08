//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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
 * Development playground for constexpression sequence
 * author: pauletti
 * date: Aug 29, 2014
 *
 */

#include "../tests.h"
#include <utility>


template <typename Type, Type ...Indices>
constexpr
auto make_index_array(std::integer_sequence<Type, Indices...>)
-> std::array<Type, sizeof...(Indices)>
{
    return std::array<Type, sizeof...(Indices)>{Indices...};
}

template<size_t N>
constexpr
auto sequence()
-> std::array<size_t, N>
{
    return make_index_array(std::make_index_sequence<N>());
}

template<int dim>
class A
{
public:
    static constexpr  std::array<size_t, dim> dims = sequence<dim>();
};

template<int dim>
constexpr  std::array<size_t, dim> A<dim>::dims;

int main()
{

    for (auto a : A<20>::dims)
        out << a << " ";
    out << endl;

    return 0;
}
