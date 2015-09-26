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
 *  Test for the l2_projection function
 *
 *  author: pauletti
 *  date: 2014-06-18
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/formula_function.h>
#include <igatools/functions/identity_function.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/space_tools.h>

#include <igatools/io/writer.h>

// TODO (pauletti, Nov 13, 2014):  add this function as the p-norm function
// in the library

template <int dim, int range, int rank>
class TestFunc : public FormulaFunction<dim, 0, range, rank>
{
private:
  using base_t   = Function<dim, 0, range, rank>;
  using parent_t = FormulaFunction<dim, 0, range, rank>;
  using self_t = TestFunc<dim, range, rank>;
  using typename base_t::GridType;
public:
  using typename parent_t::Point;
  using typename parent_t::Value;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  TestFunc(std::shared_ptr<GridType> grid)
    : parent_t(grid, IdentityFunction<dim>::const_create(grid))
  {}

  static std::shared_ptr<base_t>
  const_create(std::shared_ptr<GridType> grid)
  {
    return std::shared_ptr<base_t>(new self_t(grid));
  }

  std::shared_ptr<base_t> clone() const override
  {
    return std::make_shared<self_t>(self_t(*this));
  }

  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override
  {
    auto pt = points.begin();
    auto val = values.begin();

    for (; pt != points.end(); ++pt, ++val)
    {
      for (int i=0; i<dim; ++i)
        (*val)[i] = (*pt)[i];
      (*val)[dim] = pt->norm_square();
    }
  }


  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override
  {}

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const override
  {}
};



template<int dim, int range=1, int rank = 1, LAPack la_pack>
void test_proj(const int p, const int n_knots = 4)
{
  using Space = BSplineSpace<dim,range,rank> ;
  using RefSpace = ReferenceSpace<dim,range,rank> ;
  using Func = TestFunc<dim,range, rank>;

  auto grid = Grid<dim>::const_create(n_knots);
  auto space = Space::const_create(p, grid);

  const int n_qp = 4;
  QGauss<dim> quad(n_qp);

  auto f = Func::const_create(grid);
  auto proj_func = space_tools::projection_l2<RefSpace,la_pack>(f, space, quad);
  proj_func->print_info(out);

}


int main()
{
#if defined(USE_TRILINOS)
  const auto la_pack = LAPack::trilinos_epetra;
#elif defined(USE_PETSC)
  const auto la_pack = LAPack::petsc;
#endif

  test_proj<2, 3, 1, la_pack>(1);

  return 0;
}

