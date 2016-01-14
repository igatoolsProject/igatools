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


#include <igatools/functions/grid_function.h>
#include <igatools/functions/ig_coefficients.h>
#include <igatools/basis_functions/reference_space_basis.h>
#include <igatools/linear_algebra/epetra_vector.h>

IGA_NAMESPACE_OPEN


/**
 * @brief GridFunction built as linear combination of basis functions from ReferenceBasis
 *
 * @ingroup serializable
 */
template<int dim, int range>
class IgGridFunction :
  public GridFunction<dim, range>
{
private:
  using parent_t = GridFunction<dim, range>;
  using self_t = IgGridFunction<dim, range>;
protected:
  using typename parent_t::GridType;
  using ElementHandler = typename parent_t::ElementHandler;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  using RefBasis = ReferenceBasis<dim,range,1>;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  IgGridFunction() = default;


  virtual ~IgGridFunction() = default;

protected:
  IgGridFunction(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
                 const IgCoefficients &coeffs,
                 const std::string &dofs_property);

  IgGridFunction(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
                 const EpetraTools::Vector &coeff,
                 const std::string &dofs_property);

public:
  std::unique_ptr<ElementHandler>
  create_cache_handler() const override final;

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const RefBasis> &ref_basis,
               const IgCoefficients &coeffs,
               const std::string &dofs_property = DofProperties::active);

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<RefBasis> &ref_basis,
         const IgCoefficients &coeffs,
         const std::string &dofs_property = DofProperties::active);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const RefBasis> &ref_basis,
               const EpetraTools::Vector &coeffs,
               const std::string &dofs_property = DofProperties::active);

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<RefBasis> &ref_basis,
         const EpetraTools::Vector &coeffs,
         const std::string &dofs_property = DofProperties::active);


  virtual void print_info(LogStream &out) const override final;

  template <int sdim>
  std::shared_ptr<const IgGridFunction<sdim,range> >
  get_sub_function(const int s_id,
                   const std::shared_ptr<const Grid<sdim>> &sub_grid) const
  {
    static_assert(sdim == 0 || (sdim > 0 && sdim < dim),
                  "The dimensionality of the sub_grid is not valid.");


    typename RefBasis::template InterSpaceMap<sdim> dof_map;
    auto sub_ref_space = ref_basis_->template get_ref_sub_space<sdim>(s_id,dof_map,sub_grid);

    IgCoefficients sub_coeffs;
    const int n_sub_dofs = dof_map.size();
    for (int sub_dof = 0 ; sub_dof < n_sub_dofs ; ++ sub_dof)
      sub_coeffs[sub_dof] = coeffs_[dof_map[sub_dof]];

    auto sub_func = IgGridFunction<sdim,range>::const_create(sub_ref_space,sub_coeffs);

    return sub_func;
  }



private:
  SharedPtrConstnessHandler<RefBasis> ref_basis_;

  IgCoefficients coeffs_;

  std::string dofs_property_;

public:
  std::shared_ptr<const RefBasis> get_basis() const;

  const IgCoefficients &get_coefficients() const;

  const std::string &get_dofs_property() const;

private:
#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar);
  ///@}
#endif // SERIALIZATION

#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;


//  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &grid_function);
#endif // MESH_REFINEMENT


};

IGA_NAMESPACE_CLOSE


#ifdef SERIALIZATION

#include <igatools/functions/ig_grid_function.serial>

#endif // SERIALIZATION


#endif // __IG_GRID_FUNCTION_H

