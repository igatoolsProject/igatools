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

#ifndef __IG_GRID_FUNCTION_H
#define __IG_GRID_FUNCTION_H


#include <igatools/geometry/grid_function.h>
#include <igatools/functions/ig_coefficients.h>
#include <igatools/basis_functions/reference_space.h>
#include <igatools/linear_algebra/epetra_vector.h>

IGA_NAMESPACE_OPEN

template <int, int> class FormulaGridFunctionHandler;

/**
 * @brief GridFunction built as linear combination of basis functions from ReferenceSpace
 *
 * @ingroup serializable
 */
template<int dim, int space_dim>
class IgGridFunction :
  public GridFunction<dim, space_dim>
{
private:
  using parent_t = GridFunction<dim, space_dim>;
  using self_t = IgGridFunction<dim, space_dim>;
protected:
  using typename parent_t::GridType;
  using ElementHandler = FormulaGridFunctionHandler<dim, space_dim>;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  using IgSpace = ReferenceSpace<dim,space_dim,1>;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  IgGridFunction() = default;


  virtual ~IgGridFunction() = default;

protected:
  IgGridFunction(const SharedPtrConstnessHandler<IgSpace> &space,
                 const IgCoefficients &coeffs);

  IgGridFunction(const SharedPtrConstnessHandler<IgSpace> &space,
                 const EpetraTools::Vector &coeff);

public:
  std::unique_ptr<typename parent_t::ElementHandler>
  create_cache_handler() const override final;

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const IgSpace> &space,
               const IgCoefficients &coeffs);

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<IgSpace> &space,
         const IgCoefficients &coeffs);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const IgSpace> &space,
               const EpetraTools::Vector &coeffs);

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<IgSpace> &space,
         const EpetraTools::Vector &coeffs);


  virtual void print_info(LogStream &out) const override final;

  template <int sdim>
  std::shared_ptr<IgGridFunction<sdim,space_dim> >
  get_sub_function(const int s_id,
                   const std::shared_ptr<Grid<sdim>> &sub_grid) const
  {
    static_assert(sdim == 0 || (sdim > 0 && sdim < dim),
                  "The dimensionality of the sub_grid is not valid.");


    typename IgSpace::template InterSpaceMap<sdim> dof_map;
    auto sub_ref_space = ig_space_->template get_ref_sub_space<sdim>(s_id,dof_map,sub_grid);

    IgCoefficients sub_coeffs;
    const int n_sub_dofs = dof_map.size();
    for (int sub_dof = 0 ; sub_dof < n_sub_dofs ; ++ sub_dof)
      sub_coeffs[sub_dof] = coeffs_[dof_map[sub_dof]];

    auto sub_func = IgGridFunction<sdim,space_dim>::create(sub_ref_space,sub_coeffs);

    return sub_func;
  }



private:
  SharedPtrConstnessHandler<IgSpace> ig_space_;

  IgCoefficients coeffs_;

public:
  std::shared_ptr<const IgSpace> get_ig_space() const;

  const IgCoefficients &get_coefficients() const;

private:
#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar)
  {
    using std::to_string;
    const std::string base_name =
      "GridFunction_" + to_string(dim) + "_" +
      to_string(space_dim);

    ar &make_nvp(base_name,base_class<parent_t>(this));

    ar &make_nvp("ig_space_",ig_space_);

    ar &make_nvp("coeffs_",coeffs_);
  }
  ///@}
#endif // SERIALIZATION

};

IGA_NAMESPACE_CLOSE


#ifdef SERIALIZATION

#include <igatools/functions/ig_grid_function.serial>

#endif // SERIALIZATION


#endif // __IG_GRID_FUNCTION_H

