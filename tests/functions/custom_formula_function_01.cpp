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
 *  Test for a custom Formulafunction
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
class CustomFunction
  : public FormulaFunction<dim,codim,1,1>
{

public:
  using base_t = Function<dim, codim, 1, 1>;
  using parent_t = FormulaFunction<dim, codim, 1, 1>;
  using self_t = CustomFunction<dim, codim>;
  using typename base_t::DomainType;
  using typename parent_t::Point;
  using typename parent_t::Value;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<DomainType> &domain)
  {
    return std::shared_ptr<self_t>(new
                                   self_t(SharedPtrConstnessHandler<DomainType>(domain)));
  }

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const DomainType> &domain)
  {
    return std::shared_ptr<self_t>(new
                                   self_t(SharedPtrConstnessHandler<DomainType>(domain)));
  }

  CustomFunction(const self_t &) = default;

  virtual ~CustomFunction() = default;

  virtual void print_info(LogStream &out) const override final
  {
    Assert(false,ExcNotImplemented());
  }


protected:
  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain)
    :
    parent_t(domain)
  {}

private:
  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override
  {
    const int sp_dim = dim+codim;

    int pt = 0;
    for (auto &val : values)
    {
      const auto &point = points[pt];

      val[0] = 0.0;
      for (int i = 0 ; i < sp_dim ; ++i)
        val[0] += point[i];

      ++pt;
    }
  }

  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override
  {
    const int sp_dim = dim+codim;

//    int pt = 0;
    for (auto &val : values)
    {
//      const auto &point = points[pt];

      for (int i = 0 ; i < sp_dim ; ++i)
        val[i][0] = 1.0;

//      ++pt;
    }
  }

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const override
  {
//    Assert(false,ExcNotImplemented());
  }
};


template<int dim, int codim>
void test_custom_function(const int s_id)
{
  using std::to_string;
  out.begin_item("test_custom_function<dim=" + to_string(dim) +
                 ",codim=" + to_string(codim)+ ",range=1");

  using Function = CustomFunction<dim,codim>;


  const int sp_dim = dim+codim;

#if 0
  auto grid = Grid<dim>::create(3);

  using GridF = grid_functions::LinearGridFunction<dim,sp_dim>;
  using Grad = typename GridF::Gradient;
  using Val = typename GridF::Value;
  Grad A;
  Val b;

  const auto &sub_elem = UnitElement<sp_dim>::template get_elem<dim>(s_id);
  const auto &active_dirs = sub_elem.active_directions;
  const auto &constant_dirs = sub_elem.constant_directions;
  const auto &constant_vals = sub_elem.constant_values;
  int i = 0;
  for (const auto active_dir : active_dirs)
  {
    A[i][active_dir] = 1.0;
    ++i;
  }

  i = 0;
  for (const auto constant_dir : constant_dirs)
  {
    b[constant_dir] = constant_vals[i];
    ++i;
  }
  auto grid_func = GridF::create(grid,A,b);

  auto domain = Domain<dim,codim>::create(grid_func);
#endif

//#if 0
  auto sup_grid = Grid<sp_dim>::create(3);

  using SubGridElemMap = typename Grid<sp_dim>::template SubGridMap<dim>;
  SubGridElemMap grid_elem_map;
  LogStream myout;


  auto grid = sup_grid->template get_sub_grid<dim>(s_id,grid_elem_map);
  /*
  myout.begin_item("grid_elem_map ID: " + std::to_string(s_id));
  grid_elem_map.print_info(myout);
  myout.end_item();
  //*/

  auto sup_domain =
    Domain<sp_dim,0>::create(
      grid_functions::IdentityGridFunction<sp_dim>::create(sup_grid));
  const auto domain = sup_domain->template get_sub_domain<dim>(s_id,grid_elem_map,grid);
//#endif

  /*
    std::string filename_grid = "new_grid_s_id_" + std::to_string(s_id);
    Writer<dim,codim> writer_grid(grid,5);
    writer_grid.save(filename_grid);

    std::string filename_domain = "new_domain_s_id_" + std::to_string(s_id);
    Writer<dim,codim> writer_domain(domain,5);
    writer_domain.save(filename_domain);
  //*/
  auto F = Function::const_create(domain);

  out.begin_item("SubElem ID: " + std::to_string(s_id));
  function_values(*F);
  out.end_item();

  out.end_item();
}




int main()
{
//  test_custom_function<1, 0, 1>();
//  test_custom_function<2, 0, 2>();
//  test_custom_function<3, 0, 3>();


  for (int s_id = 0 ; s_id < 4 ; ++s_id)
  {
    test_custom_function<1,1>(s_id);
  }

  return 0;
}

