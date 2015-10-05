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


#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/bernstein_basis.h>

#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
#include <numeric>
#include <memory>

using std::reverse;
using std::accumulate;
using std::sort;

using std::shared_ptr;
using std::make_shared;



IGA_NAMESPACE_OPEN


//#define NOT_OPTIMIZED

namespace
{

#ifndef NOT_OPTIMIZED

/**
 * @class DerivativeEvaluationSymmetryManager
 *
 * This class is used to manage the symmetries in the evaluation of a k-th order
 *  derivative of a scalar function \f$ f \mathbb{R}^d -to \mathbb{R} \f$
 *  according with the Schwarz theorem (equality of mixed partials).
 *
 * The total number of derivatives is given by \f$ d^k \f$ where the number of
 * different values
 * is \f$ \binorm{d-1+k}{d-1} = \binorm{d-1+k}{k} \f$.
 *
 */
// TODO (pauletti, Nov 1, 2013): Document how this class work and its internals
template<int dim, int order>
class DerivativeSymmetryManager
{
public:
  static const int num_entries_total = pow(dim,order);
  static const int num_entries_eval = constexpr_binomial_coefficient(dim-1+order,order);
  static const int num_entries_copy = num_entries_total - num_entries_eval;

  typedef Derivatives<dim,1,1,order> Derivative_t;

  /** @name Constructors */
  ///@{
  /** Constructor */
  DerivativeSymmetryManager();

  /** Copy constructor */
  DerivativeSymmetryManager(
    const DerivativeSymmetryManager<dim,order> &in) = default;

  /** Move constructor */
  DerivativeSymmetryManager(
    DerivativeSymmetryManager<dim,order> &&in) = default;

  ~DerivativeSymmetryManager() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator */
  DerivativeSymmetryManager<dim,order> &operator=(
    const DerivativeSymmetryManager<dim,order> &) = default;

  /** Move assignment operator */
  DerivativeSymmetryManager<dim,order> &operator=(
    DerivativeSymmetryManager<dim,order> &&) = default;
  ///@}


  static const int new_deriv_order = order>0?order:1;
  const SafeSTLArray<int,num_entries_eval> &get_entries_flat_id_evaluate() const
  {
    return entries_flat_id_evaluate_;
  }

  const SafeSTLArray<int,num_entries_copy> &get_entries_flat_id_copy_to() const
  {
    return entries_flat_id_copy_to_;
  }

  const SafeSTLArray<int,num_entries_copy> &get_entries_flat_id_copy_from() const
  {
    return entries_flat_id_copy_from_;
  }

  const SafeSTLArray<TensorIndex<new_deriv_order>,num_entries_total> &get_entries_tensor_id() const
  {
    return entries_tensor_id_;
  }
private:

  bool test_if_evaluate(const SafeSTLArray<int,order> &tensor_index) const;

  /** Tensor ids of all the entries */

  SafeSTLArray<TensorIndex<new_deriv_order> ,num_entries_total> entries_tensor_id_;

  /** Flat ids of the entries that need to be computed */
  SafeSTLArray<int,num_entries_eval> entries_flat_id_evaluate_;

  /** Flat ids of the destination entries that need to be copied */
  SafeSTLArray<int,num_entries_copy> entries_flat_id_copy_to_;

  /** Flat ids of the source entries that need to be copied */
  SafeSTLArray<int,num_entries_copy> entries_flat_id_copy_from_;

  TensorSize<order> size_deriv_index_;
};




template<int dim, int order>
DerivativeSymmetryManager<dim,order>::
DerivativeSymmetryManager()
{
  size_deriv_index_.fill(dim);

  typedef MultiArrayUtils<order> MAUtils;


  Derivatives<dim,1,1,order> derivative;

  int eval_id = 0;
  int copy_id = 0;
  if (order == 0)
  {
    for (int flat_id = 0; flat_id < num_entries_total; ++flat_id)
    {
      entries_flat_id_evaluate_[eval_id] = flat_id;
      eval_id++;
    }
  }
  else
  {
    auto weights = MAUtils::compute_weight(size_deriv_index_);

    for (Index flat_id = 0; flat_id < num_entries_total; ++flat_id)
    {
      entries_tensor_id_[flat_id] = derivative.flat_to_tensor_index(flat_id);

      auto tensor_id = MAUtils::flat_to_tensor_index(flat_id,weights);

      if (test_if_evaluate(tensor_id))
      {
        entries_flat_id_evaluate_[eval_id] = flat_id;
        eval_id++;
      }
      else
      {
        entries_flat_id_copy_to_[copy_id] = flat_id;
        sort(tensor_id.begin(),tensor_id.end());
        reverse(tensor_id.begin(),tensor_id.end());
        entries_flat_id_copy_from_[copy_id] = MAUtils::tensor_to_flat_index(tensor_id,weights);
        copy_id++;
      }
    }
  }
  Assert(eval_id == num_entries_eval, ExcDimensionMismatch(eval_id,num_entries_eval));
  Assert(copy_id == num_entries_copy, ExcDimensionMismatch(copy_id,num_entries_copy));

}


template<int dim, int order>
inline
bool
DerivativeSymmetryManager<dim,order>::
test_if_evaluate(const SafeSTLArray<int,order> &tensor_index) const
{
  bool test_result = true;
  for (int i = 0; i < order-1; ++i)
  {
    if (tensor_index[i+1] > tensor_index[i])
    {
      test_result = false;
      break;
    }
  }
  return test_result;
}


template<int order>
class DerivativeSymmetryManager<0,order>
{
public:
  static const int num_entries_total = 1;
  static const int num_entries_eval  = 1;
  static const int num_entries_copy  = num_entries_total - num_entries_eval;

  typedef Derivatives<0,1,1,order> Derivative_t;

  DerivativeSymmetryManager()
  {
    entries_tensor_id_[0] = {{}};
    entries_flat_id_evaluate_[0] = 0;
  }

  const SafeSTLArray<int,num_entries_eval> &get_entries_flat_id_evaluate() const
  {
    return entries_flat_id_evaluate_;
  }

  const SafeSTLArray<int,num_entries_copy> &get_entries_flat_id_copy_to() const
  {
    return entries_flat_id_copy_to_;
  }

  const SafeSTLArray<int,num_entries_copy> &get_entries_flat_id_copy_from() const
  {
    return entries_flat_id_copy_from_;
  }

  const SafeSTLArray<TensorIndex<order>,num_entries_total> &get_entries_tensor_id() const
  {
    return entries_tensor_id_;
  }

private:
  /** Tensor ids of all the entries */
  SafeSTLArray<TensorIndex<order>,num_entries_total> entries_tensor_id_;

  /** Flat ids of the entries that need to be computed */
  SafeSTLArray<int,num_entries_eval> entries_flat_id_evaluate_;

  /** Flat ids of the destination entries that need to be copied */
  SafeSTLArray<int,num_entries_copy> entries_flat_id_copy_to_;

  /** Flat ids of the source entries that need to be copied */
  SafeSTLArray<int,num_entries_copy> entries_flat_id_copy_from_;

};


};
#endif



template <int dim, int range, int rank>
BSplineElement<dim, range, rank>::
BSplineElement(const std::shared_ptr<ContainerType> space,
               const ListIt &index,
               const PropId &prop)
  :
  parent_t(space,index,prop),
  grid_elem_(space->get_ptr_const_grid()->create_element(index,prop))
{}






template <int dim, int range, int rank>
auto
BSplineElement<dim, range, rank>::
get_bspline_space() const -> std::shared_ptr<const Space>
{
  auto bsp_space = std::dynamic_pointer_cast<const Space>(this->space_);
  Assert(bsp_space != nullptr,ExcNullPtr());
  return bsp_space;
}



template <int dim, int range, int rank>
void
BSplineElement<dim, range, rank>::
print_info(LogStream &out) const
{
  using std::to_string;

  out.begin_item("ReferenceElement<" +
                 to_string(dim) + "," +
                 to_string(range) + "," +
                 to_string(rank) + ">");
  parent_t::print_info(out);
  out.end_item();

  out.begin_item("GridElement<" + to_string(dim) + ">");
  grid_elem_->print_info(out);
  out.end_item();
}


template <int dim, int range, int rank>
void
BSplineElement<dim, range, rank>::
print_cache_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("BSplineElement<" +
                 to_string(dim) + "," +
                 to_string(range) + "," +
                 to_string(range) + "> cache:");

  out.begin_item("Splines 1D table:");
  for (int sdim = 0 ; sdim <= dim ; ++sdim)
  {
    out.begin_item("Sub-element dimension: " + to_string(sdim));
    all_splines_1D_table_[sdim].print_info(out);
    out.end_item();
  }
  out.end_item();

  out.begin_item("SpaceElement's cache:");
  SpaceElement<dim,0,range,rank,Transformation::h_grad>::print_cache_info(out);
  out.end_item();

  out.end_item();
}



IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element.inst>


