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

#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/nurbs_handler.h>
//#include <igatools/basis_functions/basis_tools.h>

//#include <igatools/base/sub_function.h>
//#include <igatools/base/exceptions.h>

#ifdef USE_NURBS


using std::to_string;
using std::endl;
using std::shared_ptr;
using std::make_shared;
using std::bind;
using std::placeholders::_1;
using std::placeholders::_2;

IGA_NAMESPACE_OPEN



template <int dim_, int range_, int rank_>
NURBS<dim_, range_, rank_>::
NURBS(const SharedPtrConstnessHandler<BSpBasis> &bsp_basis,
      const SharedPtrConstnessHandler<WeightFunction> &weight_func)
  :
  bsp_basis_(bsp_basis),
  weight_func_(weight_func)
{
#ifndef NDEBUG
  Assert(this->get_grid() == weight_func_->get_grid(),ExcMessage("Mismatching grids."));

  const auto w_as_ig_func =
    std::dynamic_pointer_cast<const IgGridFunction<dim,1>>(weight_func_.get_ptr_const_data());
  Assert(w_as_ig_func != nullptr,ExcNullPtr());

  const auto w_func_basis = w_as_ig_func->get_basis();
  Assert(w_func_basis->is_bspline(),
         ExcMessage("The basis for the weight function is not BSpline."));

  const auto &n_basis_table = this->get_spline_space()->get_dof_distribution()->get_num_dofs_table();
  int comp_id = 0;
  for (const auto &n_basis_comp : n_basis_table)
  {
    Assert(n_basis_comp == w_func_basis->get_spline_space()->get_dof_distribution()->get_num_dofs_table()[0],
           ExcMessage("Mismatching number of basis functions and weight "
                      "coefficients for scalar component " + to_string(comp_id)));

    ++comp_id;
  }
#endif
}



template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create(const std::shared_ptr<BSpBasis> &bs_basis,
       const std::shared_ptr<WeightFunction> &weight_func) -> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<BSpBasis>(bs_basis),
  SharedPtrConstnessHandler<WeightFunction>(weight_func)));
  Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif // MESH_REFINEMENT

  return sp;
}

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
const_create(const std::shared_ptr<const BSpBasis> &bs_basis,
             const std::shared_ptr<const WeightFunction> &weight_func) -> shared_ptr<const self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<BSpBasis>(bs_basis),
  SharedPtrConstnessHandler<WeightFunction>(weight_func)));
  Assert(sp != nullptr, ExcNullPtr());

  return sp;
}


template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_this_basis() const -> std::shared_ptr<const self_t >
{
  auto ref_bs = const_cast<self_t *>(this)->shared_from_this();
  auto nrb_basis = std::dynamic_pointer_cast<self_t>(ref_bs);
  Assert(nrb_basis != nullptr,ExcNullPtr());

  return nrb_basis;
}




template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_element_begin(const PropId &property) const
-> std::unique_ptr<BasisElement<dim_,0,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

//  const auto &id_elems_with_property = this->get_grid()->get_elements_with_property(property);
  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_begin(property),
  weight_func_->create_element_begin(property));
}

template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_element_end(const PropId &property) const
-> std::unique_ptr<BasisElement<dim_,0,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

//  const auto &id_elems_with_property = this->get_grid()->get_elements_with_property(property);
  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_end(property),
  weight_func_->create_element_end(property));
}


template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_ref_element_begin(const PropId &property) const
-> std::unique_ptr<ReferenceBasisElement<dim_,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_begin(property),
  weight_func_->create_element_begin(property));
}

template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_ref_element_end(const PropId &property) const
-> std::unique_ptr<ReferenceBasisElement<dim_,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_end(property),
  weight_func_->create_element_end(property));
}

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_weight_func() const -> std::shared_ptr<const WeightFunction>
{
  return weight_func_.get_ptr_const_data();
}

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_grid() const -> std::shared_ptr<const Grid<dim_>>
{
  return bsp_basis_->get_grid();
}




template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_bspline_basis() const -> const std::shared_ptr<const BSpBasis>
{
  return bsp_basis_.get_ptr_const_data();
}





template<int dim_, int range_, int rank_>
template<int sdim>
auto
NURBS<dim_, range_, rank_>::
get_sub_nurbs_basis(const int s_id,
                    InterBasisMap<sdim> &dof_map,
                    const std::shared_ptr<const Grid<sdim>> &sub_grid) const
-> std::shared_ptr<const NURBS<sdim,range_,rank_> >
{
  static_assert(sdim == 0 || (sdim > 0 && sdim < dim_),
  "The dimensionality of the sub_grid is not valid.");

  auto sub_bsp_basis = bsp_basis_->template get_sub_bspline_basis<sdim>(s_id,dof_map,sub_grid);
  auto basis_sub_grid = sub_bsp_basis->get_grid();

  auto sub_w_func = weight_func_->template get_sub_function<sdim>(s_id,basis_sub_grid);

  auto sub_nrb_basis = NURBS<sdim,range_,rank_>::const_create(sub_bsp_basis,sub_w_func);

  return sub_nrb_basis;
}


#if 0
template<int dim_, int range_, int rank_>
template<int k>
auto
NURBS<dim_, range_, rank_>::
get_sub_basis(const int s_id, InterBasisMap<k> &dof_map,
              SubGridMap<k> &elem_map) const
-> std::shared_ptr<SubBasis<k> >
{
  //TODO (martinelli Nov 27,2014): implement this function
  static_assert(k == 0 || (k > 0 && k < dim_),
  "The dimensionality of the sub_grid is not valid.");

#if 0
  using SubMap = SubMapFunction<k, dim, space_dim>;
  auto grid =  this->get_grid();
//    typename GridType::template InterGridMap<k> elem_map;
//    auto sub_grid = this->get_grid()->template get_sub_grid<k>(s_id, elem_map);

  auto sub_ref_basis = get_ref_sub_basis(s_id, dof_map, sub_grid);
  auto F = IdentityFunction<dim>::create(grid);
  auto sub_map_func = SubMap::create(sub_grid, F, s_id, *elem_map);
  auto sub_basis = SubBasis<k>::create(sub_ref_basis, sub_map_func);
  return sub_basis;
#endif
  Assert(false,ExcNotImplemented());
  AssertThrow(false,ExcNotImplemented());
  return nullptr;
}
#endif

template <int dim_, int range_, int rank_>
void
NURBS<dim_, range_, rank_>::
print_info(LogStream &out) const
{
  out.begin_item("BSpline Basis:");
  bsp_basis_->print_info(out);
  out.end_item();

  out.begin_item("Weight function :");
  weight_func_->print_info(out);
  out.end_item();
}



template <int dim_, int range_, int rank_>
bool
NURBS<dim_, range_, rank_>::
is_bspline() const
{
  return false;
}

#if 0
template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_degree_table() const -> const DegreeTable &
{
  return this->bsp_basis_->get_degree_table();
}
#endif





#if 0
template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_periodicity() const -> const PeriodicityTable &
{
  return bsp_basis_->get_periodicity();
}
#endif

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_end_behaviour_table() const -> const EndBehaviourTable &
{
  return bsp_basis_->get_end_behaviour_table();
}


template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_cache_handler() const -> std::unique_ptr<BasisHandler<dim_,0,range_,rank_>>
{
  return std::unique_ptr<Handler>(new Handler(this->get_this_basis()));
}


template<int dim_, int range_, int rank_>
std::shared_ptr<const SplineSpace<dim_,range_,rank_>>
                                                   NURBS<dim_, range_, rank_>::
                                                   get_spline_space() const
{
  return bsp_basis_->get_spline_space();
}


#ifdef MESH_REFINEMENT


template<int dim_, int range_, int rank_>
void
NURBS<dim_, range_, rank_>::
refine_h(const Size n_subdivisions)
{
  //the refinement of the BSpline also refines the weight_fucntion (they share the same Grid)
  bsp_basis_.get_ptr_data()->refine_h(n_subdivisions);
}

template<int dim_, int range_, int rank_>
void
NURBS<dim_, range_, rank_>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert,
  const Grid<dim_> &old_grid)
{
  auto bsp_basis_pre_refinement =
    std::dynamic_pointer_cast<const BSpBasis>(bsp_basis_->get_basis_previous_refinement());
  Assert(bsp_basis_pre_refinement != nullptr,ExcNullPtr());

  auto weight_func_pre_refinement_ =
    std::dynamic_pointer_cast<const WeightFunction>(weight_func_->get_grid_function_previous_refinement());
  Assert(weight_func_pre_refinement_ != nullptr,ExcNullPtr());


  this->ref_basis_previous_refinement_ =
    self_t::const_create(bsp_basis_pre_refinement,weight_func_pre_refinement_);
}




#endif //MESH_REFINEMENT


#ifdef SERIALIZATION

template<int dim_, int range_, int rank_>
template<class Archive>
void
NURBS<dim_, range_, rank_>::
serialize(Archive &ar)
{
  using std::to_string;
  const std::string base_name = "ReferenceBasis_" +
                                to_string(dim_) + "_" +
                                to_string(0) + "_" +
                                to_string(range_) + "_" +
                                to_string(rank_);

  ar &make_nvp(base_name,base_class<RefBasis>(this));
  ar &make_nvp("bsp_basis_",bsp_basis_);

  ar &make_nvp("weight_func_",weight_func_);
}
///@}

#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs.inst>

#endif // #ifdef USE_NURBS

