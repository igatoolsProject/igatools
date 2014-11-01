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
 *  Test for the BSplineSpace element iterator to get the face values
 *  author: pauletti
 *  date: 2013-10-23
 */

#include "../tests.h"

#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/quadrature_lib.h>


template<int dim, int k=dim-1>
void sub_elem_values(const int n_knots, const int deg)
{
	OUTSTART

	auto grid = CartesianGrid<dim>::create(n_knots);
	using Space = NewBSplineSpace<dim>;
	using ElementHandler = typename Space::ElementHandler;
	auto space = Space::create(deg, grid);

	const int n_qp = 2;
	QGauss<k>   k_quad(n_qp);
	QGauss<dim> quad(n_qp);
	auto flag = NewValueFlags::value;//|NewValueFlags::gradient|NewValueFlags::hessian;
	ElementHandler cache(space);
	cache.template reset<k>(flag, k_quad);
	cache.template reset<dim>(flag, quad);
	auto elem = space->begin();
	auto end =  space->end();
	cache.template init_cache<k>(elem);
	cache.template init_cache<dim>(elem);
	for (; elem != end; ++elem)
	{
		if (elem->is_boundary())
		{
			cache.template fill_cache<dim>(elem, 0);
			out << "Element" << elem->get_flat_index() << endl;
			elem->template get_values<0,dim>(0).print_info(out);
			for(auto &s_id : UnitElement<dim>::template elems_ids<k>())
			{
				if (elem->is_boundary(s_id))
				{
					cache.template fill_cache<k>(elem, s_id);
					out << "Sub Element: " << s_id << endl;
					out.begin_item("Values basis functions:");
					auto values = elem->template get_values<0,k>(s_id);
					values.print_info(out);
					out.end_item();
				}
			}
		}
	}
	OUTEND
}
//	auto values    = elem->template get_values<0,k>(0);
//	auto gradients = elem->template get_values<1,k>(0);
//	auto hessians  = elem->template get_values<2,k>(0);

//	out.begin_item("Values basis functions:");
//	values.print_info(out);
//	out.end_item();

//	out.begin_item("Gradients basis functions:");
//	gradients.print_info(out);
//	out.end_item();
//
//	out.begin_item("Hessians basis functions:");
//	hessians.print_info(out);
//	out.end_item();


//    for (; elem != end; ++elem)
//    {
//        if (elem->is_boundary())
//        {
//            elem->fill_cache();
//            out << "Element" << elem->get_flat_index() << endl;
//            elem->get_basis_values().print_info(out);
//
//            for (int face = 0; face < num_faces; ++face)
//            {
//                if (elem->is_boundary(face))
//                {
//                    elem->fill_face_cache(face);
//                    out << "Face " << face << endl;
//                    out << "values: " << endl;
//                    elem->get_face_basis_values(face).print_info(out);
//                    out << endl;
//                }
//            }
//        }
//    }
//    out << endl;
//}




int main()
{
    out.depth_console(10);
    sub_elem_values<2,1>(2,1);

    return 0;
}

