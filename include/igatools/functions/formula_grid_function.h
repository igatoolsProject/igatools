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

#ifndef __FORMULA_GRID_FUNCTION_H_
#define __FORMULA_GRID_FUNCTION_H_

#include <igatools/base/value_types.h>
#include <igatools/functions/grid_function.h>
#include <igatools/functions/grid_function_element.h>
#include <igatools/functions/sub_grid_function.h>

IGA_NAMESPACE_OPEN

template <int, int> class FormulaGridFunctionHandler;

/**
 *
 */
template<int dim, int range>
class FormulaGridFunction :
  public GridFunction<dim, range>
{
private:
  using parent_t =  GridFunction<dim, range>;
  using self_t = FormulaGridFunction<dim, range>;
protected:
  using typename parent_t::GridType;
  using Handler = FormulaGridFunctionHandler<dim, range>;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  FormulaGridFunction(const SharedPtrConstnessHandler<GridType> &grid);

  virtual ~FormulaGridFunction() = default;

  std::unique_ptr<typename parent_t::Handler>
  create_cache_handler() const;


  template <int sdim>
  using SubGridElemMap = typename Grid<dim>::template SubGridMap<sdim>;

  template <int sdim>
  std::shared_ptr<const SubGridFunction<sdim,dim,range> >
  get_sub_function(const int s_id,
                   const SubGridElemMap<sdim> &sub_grid_elem_map,
                   const std::shared_ptr<const Grid<sdim>> &sub_grid) const
  {
    static_assert(sdim == 0 || (sdim > 0 && sdim < dim),
                  "The dimensionality of the sub_grid is not valid.");

    auto sub_func = SubGridFunction<sdim,dim,range>::const_create(
                      this->shared_from_this(),s_id,sub_grid_elem_map,sub_grid);
//    AssertThrow(false,ExcNotImplemented());
    /*
        typename RefBasis::template InterSpaceMap<sdim> dof_map;
        auto sub_ref_space = ref_basis_->template get_ref_sub_space<sdim>(s_id,dof_map,sub_grid);

        IgCoefficients sub_coeffs;
        const int n_sub_dofs = dof_map.size();
        for (int sub_dof = 0 ; sub_dof < n_sub_dofs ; ++ sub_dof)
          sub_coeffs[sub_dof] = coeffs_[dof_map[sub_dof]];

        auto sub_func = IgGridFunction<sdim,range>::const_create(sub_ref_space,sub_coeffs);

        return sub_func;
      //*/
    return sub_func;
  }

public:

  virtual void evaluate_0(const ValueVector<GridPoint> &points,
                          ValueVector<Value> &values) const = 0;

  virtual void evaluate_1(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<1>> &values) const = 0;

  virtual void evaluate_2(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<2>> &values) const = 0;


private:



};

IGA_NAMESPACE_CLOSE

#endif
