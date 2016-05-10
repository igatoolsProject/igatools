//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
 *  Test for creation of basis only passing a domain to a function
 *
 *  author: pauletti
 *  date: 2016/05/02
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/basis_functions/physical_basis.h>



template <int dim, int range=1, int rank=1, int codim = 0>
std::shared_ptr<const PhysicalBasis<dim,range,rank,codim> >
create_basis(std::shared_ptr<Domain<dim,codim> > domain)
{
  OUTSTART

  auto grid = domain->get_grid_function()->get_grid();
  const int deg = 2;
  auto ref_basis = BSpline<dim,range,rank>::create(
                     SplineSpace<dim,range,rank>::create(deg,grid));



  using PhysBasis = PhysicalBasis<dim,range,rank,codim>;
  auto phys_basis = PhysBasis::const_create(ref_basis,domain,Transformation::h_grad);


  OUTEND

  return phys_basis;
}

//template <int dim, int range=1, int rank=1, int codim = 0>
//std::shared_ptr<const PhysicalBasis<dim,range,rank,codim> >
//create_const_basis(std::shared_ptr<const Domain<dim,codim> > domain)
//{
//  OUTSTART
//
//  auto grid = domain->get_grid_function()->get_grid();
//  const int deg = 2;
//  auto ref_basis = BSpline<dim,range,rank>::const_create(
//                     SplineSpace<dim,range,rank>::const_create(deg,grid));
//
//
//
//  using PhysBasis = PhysicalBasis<dim,range,rank,codim>;
//  auto phys_basis = PhysBasis::const_create(ref_basis,domain,Transformation::h_grad);
//
//
//  OUTEND
//
//  return phys_basis;
//}


template <int dim, int range=1, int rank=1, int codim = 0>
void
basis_through_domain()
{
	const int n_knots=4;
	auto grid  = Grid<dim>::create(n_knots);
	using GridFunc = grid_functions::BallGridFunction<dim>;
	auto grid_func = GridFunc::create(grid);

	using Domain = Domain<dim,codim>;
	auto domain = Domain::create(grid_func);

	auto phys_basis = create_phys_basis<dim,range,rank,codim>(domain);
	phys_basis->print_info(out);
}



int main()
{

  basis_through_domain<1>();

  return  0;
}
