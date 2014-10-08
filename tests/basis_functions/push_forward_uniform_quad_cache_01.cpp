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
 *  Test for developing the  PhysicalUniformQuadCache
 *
 *  author: pauletti
 *  date: Aug 21, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/push_forward_uniform_quad_cache.h>
#include <igatools/geometry/push_forward_element_accessor.h>
#include <igatools/geometry/identity_mapping.h>

template <int dim, int range=1, int rank=1>
void uniform_pf_cache(const ValueFlags flag,
                      const int n_knots = 5)
{
    OUTSTART
    using PF = PushForward<Transformation::h_grad,dim,0>;

    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto map = IdentityMapping<dim>::create(grid);
    auto push_forward = PF::create(map);

    auto quad = QGauss<dim>(2);
    PushFowardUniformQuadCache<PF> cache(push_forward, flag, quad);
    cache.print_info(out);

    OUTEND
}


template <int dim, int range=1, int rank=1>
void pf_cache_init_elem(const ValueFlags flag,
                        const int n_knots = 5)
{
    OUTSTART

    using PF = PushForward<Transformation::h_grad,dim,0>;
    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto map = IdentityMapping<dim>::create(grid);
    auto push_forward = PF::create(map);

    auto quad = QGauss<dim>(2);
    PushFowardUniformQuadCache<PF> cache(push_forward, flag, quad);

    auto elem = push_forward->begin();
    cache.init_element_cache(elem);
    elem->print_cache_info(out);

    OUTEND
}



template <int dim, int range=1, int rank=1>
void pf_cache_fill_elem(const ValueFlags flag,
                        const int n_knots = 5)
{
    OUTSTART

    using PF = PushForward<Transformation::h_grad,dim,0>;
    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto map = IdentityMapping<dim>::create(grid);
    auto push_forward = PF::create(map);

    auto quad = QGauss<dim>(2);
    PushFowardUniformQuadCache<PF> cache(push_forward, flag, quad);

    auto elem = push_forward->begin();
    auto end = push_forward->end();

    cache.init_element_cache(elem);

    for (; elem != end; ++elem)
    {
        cache.fill_element_cache(elem);
        elem->print_cache_info(out);
    }

    OUTEND
}




int main()
{
    out.depth_console(10);

    uniform_pf_cache<1>(ValueFlags::value);
    uniform_pf_cache<2>(ValueFlags::value);

    pf_cache_init_elem<1>(ValueFlags::tran_value);
    pf_cache_fill_elem<1>(ValueFlags::tran_value);

    return  0;
}
