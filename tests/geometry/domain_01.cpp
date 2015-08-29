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

/**
 *  @file
 *  @brief Developing domain class
 *  @author pauletti
 *  @date 2015
 */

#include "../tests.h"

#include <igatools/geometry/physical_domain.h>



template<int dim, int codim>
void domain()
{
    OUTSTART

	using Grid = CartesianGrid<dim>;
    using Function = Function<dim,0,dim+codim,1>;
    using Domain   = PhysicalDomain<dim, codim>;

    auto grid = Grid::const_create();
    //auto F = Function::create();
    std::shared_ptr<const Function> F;

    Domain dom(grid, F);

    OUTEND
}


int main()
{
    domain<1,0>();
    domain<2,0>();
    domain<3,0>();

    return 0;
}

