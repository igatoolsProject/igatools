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
 *  Test for BSplineBasis using a non range-homogeneous space
 *  Evaluates values gradients and derivatives at two quad point
 *  on each element without the use of the cache.
 *  This test computes the same quantities of bspline_element_iterator_04.cpp
 *
 *  author: martinelli
 *  date: May 08, 2014
 *
 */
// TODO (pauletti, Jun 2, 2014): no need to print info of the spaces
// TODO (pauletti, Jun 2, 2014): rename this file to have 04 at the end
#include "../tests.h"


#include <igatools/basis_functions/bspline.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_element.h>

using std::shared_ptr;

template <int dim>
shared_ptr<BSpline<dim,dim,1> >
create_basis(const int num_knots) ;

template <>
shared_ptr<BSpline<2,2,1> >
create_basis<2>(const int num_knots)
{
  auto knots = Grid<2>::create(num_knots);

  using Space = SplineSpace<2,2,1>;
  using Basis = BSpline<2,2,1>;
  typename Space::DegreeTable degree_table = {{{3,2}},{{2,3}}} ;

  using PeriodicityTable = typename Space::PeriodicityTable;
  using EndBehaviourTable = typename Basis::EndBehaviourTable;

  return Basis::create(
           Space::create(degree_table, knots,
                         Space::get_multiplicity_from_regularity(
                           InteriorReg::maximum,
                           degree_table,
                           knots->get_num_intervals()),
                         PeriodicityTable(true,SafeSTLArray<bool, 2>(false))),
           EndBehaviourTable(true,SafeSTLArray<BasisEndBehaviour,2>(BasisEndBehaviour::interpolatory))) ;
}

template <>
shared_ptr<BSpline<3,3,1> >
create_basis<3>(const int num_knots)
{
  auto knots = Grid<3>::create(num_knots);

  using Space = SplineSpace<3,3,1>;
  using Basis = BSpline<3,3,1>;
  typename Space::DegreeTable degree_table = { {{3,2,2}},{{2,3,2}},{{2,2,3}}} ;

  using PeriodicityTable = typename Space::PeriodicityTable;
  using EndBehaviourTable = typename Basis::EndBehaviourTable;

  return Basis::create(
           Space::create(degree_table, knots,
                         Space::get_multiplicity_from_regularity(
                           InteriorReg::maximum,
                           degree_table,
                           knots->get_num_intervals()),
                         PeriodicityTable(true,SafeSTLArray<bool,3>(false))),
           EndBehaviourTable(true,SafeSTLArray<BasisEndBehaviour,3>(BasisEndBehaviour::interpolatory))) ;
}


template< int dim>
void do_test()
{
  OUTSTART
//    out << "domain, range and rank: " << dim_domain << "," << dim_domain << ",1" << endl ;

  const int num_knots = 3 ;

  auto basis = create_basis<dim>(num_knots) ;

  const int n_points = 2;
  auto quad = Quadrature<dim>::create(QGauss<dim>(n_points).get_points()) ;

  using std::to_string;

  using Elem = typename BSpline<dim,dim,1>::ElementAccessor;
  using _Value = typename Elem::_Value;
  using _Gradient = typename Elem::_Gradient;
  using _Hessian = typename Elem::_Hessian;
  using _Divergence = typename Elem::_Divergence;


  {
    auto elem = basis->begin();
    for (; elem != basis->end(); ++elem)
    {
      out << "Element " << elem->get_index() << std::endl;
      const auto values = elem->template evaluate_basis_at_points<_Value>(quad,DofProperties::active);
      out.begin_item("Basis function values:");
      values.print_info(out);
      out.end_item();
    }
  }

  {
    auto elem = basis->begin();
    for (; elem != basis->end(); ++elem)
    {
      auto elem = basis->begin();
      for (; elem != basis->end(); ++elem)
      {
        out << "Element " << elem->get_index() << std::endl;
        const auto gradients = elem->template evaluate_basis_at_points<_Gradient>(quad,DofProperties::active);
        out.begin_item("Basis function gradients:");
        gradients.print_info(out);
        out.end_item();
      }
    }
  }

  {
    auto elem = basis->begin();
    for (; elem != basis->end(); ++elem)
    {
      auto elem = basis->begin();
      for (; elem != basis->end(); ++elem)
      {
        out << "Element " << elem->get_index() << std::endl;
        const auto hessians = elem->template evaluate_basis_at_points<_Hessian>(quad,DofProperties::active);
        out.begin_item("Basis function hessians:");
        hessians.print_info(out);
        out.end_item();
      }
    }
  }

  OUTEND

}


int main()
{
  out.depth_console(10); //to be removed after test finished

  do_test<2>() ;
  do_test<3>() ;

  return 0;
}
