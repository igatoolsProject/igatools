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

#ifndef __IG_FUNCTION_H
#define __IG_FUNCTION_H

#include <igatools/base/value_types.h>
#include <igatools/functions/function.h>
#include <igatools/functions/ig_coefficients.h>

#include <igatools/basis_functions/spline_space.h>
#include <igatools/linear_algebra/epetra_vector.h>
//#include <igatools/basis_functions/bspline.h>
//#include <igatools/basis_functions/nurbs.h>


#include <boost/fusion/include/filter_if.hpp>
//#include <boost/fusion/include/iterator.hpp>
#include <boost/fusion/include/tag_of.hpp>
#include <boost/fusion/include/key_of.hpp>
#include <boost/mpl/not_equal_to.hpp>
#include <boost/fusion/include/begin.hpp>
IGA_NAMESPACE_OPEN



template <int,int,int,int>
class PhysicalBasis;


template <int,int,int,int>
class BasisHandler;

template <int,int,int>
class BSplineHandler;

template <int,int,int>
class NURBSHandler;

template <int,int,int,int>
class PhysicalBasisHandler;

template <int,int,int,int>
class BasisElement;

template <int,int,int>
class BSplineElement;

template <int,int,int>
class NURBSElement;

template <int,int,int,int>
class PhysicalBasisElement;


/**
 * @brief Function built as linear combination of basis functions from PhysicalBasis
 *
 * @ingroup serializable
 */
template<int dim,int codim,int range,int rank>
class IgFunction :
  public Function<dim,codim,range,rank>
{
public:

  using CoeffType = IgCoefficients;


private:
  using base_t = Function<dim,codim,range,rank>;
  using parent_t = Function<dim,codim,range,rank>;
  using self_t = IgFunction<dim,codim,range,rank>;
  using PhysBasis = PhysicalBasis<dim,range,rank,codim>;

public:
  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  IgFunction() = default;

#ifdef IGATOOLS_USES_TRILINOS
  //TODO (pauletti, Mar 23, 2015): should we make this private?
  IgFunction(const SharedPtrConstnessHandler<PhysBasis> &basis,
             const EpetraTools::Vector &coeff,
             const std::string &dofs_property);
#endif // IGATOOLS_USES_TRILINOS

  IgFunction(const SharedPtrConstnessHandler<PhysBasis> &basis,
             const IgCoefficients &coeff,
             const std::string &dofs_property);



  virtual ~IgFunction() = default;

//  using typename parent_t::topology_variant;
//  using typename parent_t::eval_pts_variant;
  using typename parent_t::Point;
  using typename parent_t::Value;
  using typename parent_t::Gradient;
  using typename parent_t::ElementIterator;
  using typename parent_t::ElementAccessor;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  using typename parent_t::DomainType;

public:


  std::unique_ptr<typename parent_t::Handler>
  create_cache_handler() const override final;

#ifdef IGATOOLS_USES_TRILINOS
  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const PhysBasis> &basis,
               const EpetraTools::Vector &coeff,
               const std::string &dofs_property = DofProperties::active);

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<PhysBasis> &basis,
         const EpetraTools::Vector &coeff,
         const std::string &dofs_property = DofProperties::active);
#endif // IGATOOLS_USES_TRILINOS

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const PhysBasis> &basis,
               const IgCoefficients &coeff,
               const std::string &dofs_property = DofProperties::active);


  static std::shared_ptr<self_t>
  create(const std::shared_ptr<PhysBasis> &basis,
         const IgCoefficients &coeff,
         const std::string &dofs_property = DofProperties::active);



  std::shared_ptr<const PhysBasis> get_basis() const;

  const CoeffType &get_coefficients() const;

  const std::string &get_dofs_property() const;

  self_t &operator +=(const self_t &fun);

  virtual void print_info(LogStream &out) const override final;


  template <int sdim>
  using SubFunc = typename base_t::template SubFunc<sdim>;


  std::shared_ptr<const SubFunc<(dim>0)?dim-1:0> >
  get_sub_function(const int s_id,
                   const std::shared_ptr<const Grid<(dim>0)?dim-1:0>> &sub_grid) const override final
  {
    return get_sub_function_impl<(dim>0)?dim-1:0>(s_id,sub_grid);
  }


  template <int sdim>
  std::shared_ptr<const Function<sdim,codim+(dim-sdim),range,rank> >
  get_sub_function_impl(const int s_id,
                        const std::shared_ptr<const Grid<sdim>> &sub_grid,
                        EnableIf<((dim > 0) &&(sdim >=0))> * = nullptr) const
  {
    static_assert(sdim == 0 || (sdim > 0 && sdim < dim),
                  "The dimensionality of the sub_grid is not valid.");



    typename PhysBasis::template InterBasisMap<sdim> dof_map;
    typename Grid<dim>::template SubGridMap<sdim> elem_map;
    auto sub_basis = basis_->template get_sub_basis<sdim>(s_id,dof_map,sub_grid,elem_map);


    IgCoefficients sub_coeffs;
    const int n_sub_dofs = dof_map.size();
    for (int sub_dof = 0 ; sub_dof < n_sub_dofs ; ++ sub_dof)
      sub_coeffs[sub_dof] = coeffs_[dof_map[sub_dof]];

    auto sub_func = IgFunction<sdim,codim+(dim-sdim),range,rank>::const_create(sub_basis,sub_coeffs);

    return sub_func;
  }

  template <int sdim>
  std::shared_ptr<const Function<sdim,codim+(dim-sdim),range,rank> >
  get_sub_function_impl(const int s_id,
                        const std::shared_ptr<const Grid<sdim>> &sub_grid,
                        EnableIf<!((dim > 0) &&(sdim >=0))> * = nullptr) const
  {
    AssertThrow(false,ExcNotImplemented());
    return nullptr;
  }

private:

  SharedPtrConstnessHandler<PhysBasis> basis_;

  CoeffType coeffs_;

  std::string dofs_property_;

private:

#ifdef IGATOOLS_WITH_MESH_REFINEMENT

//  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &ig_function);

  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;

#endif // IGATOOLS_WITH_MESH_REFINEMENT

#ifdef IGATOOLS_WITH_SERIALIZATION
  /**
   * @name Functions needed for the serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar);
  ///@}
#endif // IGATOOLS_WITH_SERIALIZATION
};



IGA_NAMESPACE_CLOSE

#ifdef IGATOOLS_WITH_SERIALIZATION

#include <igatools/functions/ig_function.serial>

#endif // IGATOOLS_WITH_SERIALIZATION

#endif


