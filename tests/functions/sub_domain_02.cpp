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
 *  Test for the creation of the subdomain from a domain
 *  built using an BallGridFunction
 *  author: martinelli
 *  date: Jan 12, 2016
 */

#include "../tests.h"

#include <igatools/functions/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/io/writer.h>


#include "function_test.h"




template<int dim, int codim>
void test_sub_domain(const int s_id)
{
  using std::to_string;
  out.begin_item("test_sub_domain<dim=" + to_string(dim) +
                 ",codim=" + to_string(codim)+ ">");


  const int sp_dim = dim+codim;


  std::shared_ptr<Grid<dim>> grid;

  std::shared_ptr<const Domain<dim,codim>> domain;


  auto sup_grid = Grid<sp_dim>::create(3);

  using SubGridElemMap = typename Grid<sp_dim>::template SubGridMap<dim>;
  SubGridElemMap grid_elem_map;


  grid = sup_grid->template get_sub_grid<dim>(s_id,grid_elem_map);

  auto sup_domain =
    Domain<sp_dim,0>::create(
      grid_functions::BallGridFunction<sp_dim>::create(sup_grid));
  domain = sup_domain->template get_sub_domain<dim>(s_id,grid_elem_map,grid);


  string msg = "_dim_" + std::to_string(dim) +
               "_codim_" + std::to_string(codim) +
               "_s_id_" + std::to_string(s_id);

  /*
  std::string filename_grid = "grid" + msg;
  Writer<dim,codim> writer_grid(grid,5);
  writer_grid.save(filename_grid);
  //*/

  std::string filename_domain = "domain" + msg;
  Writer<dim,codim> writer_domain(domain,10);
  writer_domain.save(filename_domain);
  //*/
  auto F = CustomFunction<dim,codim>::const_create(domain);

  out.begin_item("SubElem ID: " + std::to_string(s_id));
  function_values(*F);

  out << "Measure = " << compute_measure<dim,codim>(*domain);

  out.end_item();

  out.end_item();
}


template <int dim>
void do_test()
{
  for (int s_id = 0 ; s_id < UnitElement<dim+1>::n_faces ; ++s_id)
  {
    test_sub_domain<dim,1>(s_id);
  }
}


int main()
{

  do_test<1>();
  do_test<2>();

  return 0;
}

