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

/* test_utils.h
 *
 *  Created on: Aug 20, 2015
 *      Author: pauletti
 */

#ifndef __TESTS_TEST_UTILS_H_
template <int dim>
auto non_uniform_grid(const int max_nodes = 5)
{
	SafeSTLArray<SafeSTLVector<Real>,dim> knots;
    for(int i=0; i<dim; ++i)
    {
    	knots[i].resize(max_nodes-i);
    	auto prog=knots[i];
    	 std::iota (prog.begin(), prog.end(),0);
    	std::partial_sum (prog.begin(), prog.end(), knots[i].begin());
    }

    return CartesianGrid<dim>::create(knots);
}

#define __TESTS_TEST_UTILS_H_





#endif /* TESTS_TEST_UTILS_H_ */
