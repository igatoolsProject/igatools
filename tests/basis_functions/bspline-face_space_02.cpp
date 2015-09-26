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
 *  Test for the BsplineSpace class subspace extraction
 *
 *  author: pauletti
 *  date: 2014-11-18
 */

#include "../tests.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>


template<int sub_dim, int dim, int range=1, int rank=1>
void sub_space(TensorSize<dim> n, const int degree = 1)
{
  OUTSTART

  using Space = BSplineSpace<dim, range, rank>;

  auto grid = Grid<dim>::const_create(n);
  auto space = Space::const_create(degree, grid);

  typename Space::template InterSpaceMap<sub_dim> dof_map;

  for (auto i : UnitElement<dim>::template elems_ids<sub_dim>())
  {
    std::map<Index,Index> elem_map;
    auto sub_grid = space->get_ptr_const_grid()->template get_sub_grid<sub_dim>(i, elem_map);
    out.begin_item(to_string(i) + "-th " + "sub space:");
    auto sub_space =
      space->template get_sub_space<sub_dim>(i, dof_map, sub_grid, elem_map);
    out.begin_item("Space:");
    sub_space->print_info(out);
    out.end_item();

    out.begin_item("Dofs sub element to space mapping:");
    dof_map.print_info(out);
    out.end_item();
    out.end_item();
  }

  OUTEND
}



int main()
{

  sub_space<0,1>(TensorSize<1>(sequence<1>(2)));
  sub_space<1,2>(TensorSize<2>(sequence<2>(2)));
  sub_space<2,3>(TensorSize<3>(sequence<3>(2)));

  return  0;
}


//#include "../tests.h"
//#include <igatools/basis_functions/bspline_space.h>
//#include <igatools/basis_functions/space_tools.h>
//
//template< int dim_domain, int dim_range, int rank >
//void run_test()
//{
//    out << "-------------------------" <<endl;
//    out << "dim_dimain = " << dim_domain << endl;
//    out << "dim_range  = " << dim_range << endl;
//    out << "rank       = " << rank << endl;
//
//    const int degree=1;
//
//    auto grid = Grid<dim_domain>::const_create(3);
//    auto space = BSplineSpace<dim_domain, dim_range, rank>::const_create(degree, grid);
//
//    SafeSTLVector<Index> dof_map;
//
//    int face_id = 0;
//
//    auto face_space = space->get_face_space(face_id,dof_map);
//
//    for (SafeSTLVector<Index>::iterator it=dof_map.begin() ; it < dof_map.end(); ++it)
//        out<< "face_id = "<< face_id << ", dof_id = "<< *it << endl;
//
//    out << "-------------------------" <<endl;
//}
//
//
//
//int main()
//{
//    //out.depth_console(10);
//
//    run_test<2,1,1>();
//    run_test<2,2,1>();
//    run_test<3,1,1>();
//
//    return  0;
//}
