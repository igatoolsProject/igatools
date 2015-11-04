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

#include <igatools/basis_functions/phys_space_element_handler.h>
#include <igatools/basis_functions/physical_space_element.h>

#include <functional>

using std::shared_ptr;



IGA_NAMESPACE_OPEN


namespace
{
#if 0
auto
pushforward_to_mapping_flag(
  const Transformation type,
  const ValueFlags flags,
  const TransformationFlags transf_flags)
-> ValueFlags
{
  ValueFlags transfer_flag =
  ValueFlags::measure |
  ValueFlags::w_measure |
  ValueFlags::outer_normal|
  ValueFlags::boundary_normal|
  ValueFlags::point|
//    ValueFlags::value|
  ValueFlags::gradient|
  ValueFlags::hessian;

  ValueFlags map_flag = flags & transfer_flag;

  if (type == Transformation::h_grad)
  {
    if (contains(transf_flags, TransformationFlags::tran_value))
    {}

    if (contains(transf_flags, TransformationFlags::tran_gradient))
      map_flag|= (ValueFlags::inv_gradient);

    if (contains(transf_flags, TransformationFlags::tran_hessian))
    {
      map_flag |= (ValueFlags::hessian | ValueFlags::inv_gradient);
    }
  }
//    else if (type == Transformation::h_div)
//    {
//        if (contains(v_flag,ValueFlags::tran_value))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_face_gradient);
//        if (contains(v_flag,ValueFlags::tran_gradient))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_hessian |
//                          ValueFlags::map_face_gradient |
//                          ValueFlags::map_face_hessian);
//        if (contains(v_flag,ValueFlags::tran_hessian))
//            AssertThrow(false,ExcNotImplemented());
//    }
//    else if (type == Transformation::h_curl)
//    {
//        AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_value))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_face_gradient);
//        if (contains(v_flag,ValueFlags::tran_gradient))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_hessian |
//                          ValueFlags::map_face_gradient |
//                          ValueFlags::map_face_hessian);
//        if (contains(v_flag,ValueFlags::tran_hessian))
//            AssertThrow(false,ExcNotImplemented());
//    }
//    else if (type == Transformation::l_2)
//    {
//        AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_value))
//            AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_gradient))
//            AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_hessian))
//            AssertThrow(false,ExcNotImplemented());
//    }
//
//
//
//    // We fill extra stuff as the computation is performed anyways
//    if (contains(fill_flag , ValueFlags::measure))
//        fill_flag |= (ValueFlags::map_gradient |
//                      ValueFlags::map_face_gradient);
//
//    if (contains(fill_flag , ValueFlags::map_inv_gradient))
//        fill_flag |= (ValueFlags::map_gradient |
//                      ValueFlags::measure |
//                      ValueFlags::map_face_gradient |
//                      ValueFlags::face_measure);
//
//    if (contains(fill_flag , ValueFlags::map_inv_hessian))
//        fill_flag |= (ValueFlags::map_hessian |
//                      ValueFlags::map_face_hessian);

  return map_flag;
}
#endif

};


template<int dim, int codim>
using MapFunc= Function<dim, 0, dim + codim, 1>;


template<int dim_,int range_,int rank_,int codim_>
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
PhysSpaceElementHandler(std::shared_ptr<const PhysSpace> space)
  :
  base_t(space),
  ref_space_handler_(space->get_reference_basis()->create_cache_handler()),
  phys_domain_handler_(space->get_physical_domain()->create_cache_handler()),
  phys_space_(space)
{}

#if 0
template<int dim_,int range_,int rank_,int codim_>
auto
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
create(std::shared_ptr<const PhysSpace> space) -> std::unique_ptr<self_t>
{
  auto handler = std::unique_ptr<self_t>(new self_t(space));
  Assert(handler != nullptr,ExcNullPtr());
  return handler;
}
#endif


#if 0
template<int dim_,int range_,int rank_,int codim_>
template<int sub_elem_dim>
void
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
ResetDispatcher::
operator()(const Quadrature<sub_elem_dim> &quad)
{
  flags_[sub_elem_dim] = flag_in_;
  ref_space_handler_.reset_selected_elements(
    space_to_ref_flag(PhysSpace::PushForwardElem::type, flag_in_),
    quad,
    elements_flat_id_);


  ValueFlags map_flags;
  TransformationFlags transf_flags;
  space_to_pf_flag(flag_in_,map_flags, transf_flags);

  auto mapping_flags = pushforward_to_mapping_flag(
                         PhysSpace::PushForwardElem::type,
                         flag_in_,
                         transf_flags);
  phys_domain_handler_.template reset<sub_elem_dim>(mapping_flags,quad);
}

template<int dim_,int range_,int rank_,int codim_>
void
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
reset_selected_elements(
  const ValueFlags &flag,
  const eval_pts_variant &eval_points,
  const SafeSTLVector<int> &elements_flat_id)
{
  auto reset_selected_elems_dispatcher =
    ResetDispatcher(flag,elements_flat_id,*ref_space_handler_,phys_domain_handler_,flags_);
  boost::apply_visitor(reset_selected_elems_dispatcher,eval_points);
}





template<int dim_,int range_,int rank_,int codim_>
template<int sub_elem_dim>
void
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
InitCacheDispatcher::
operator()(const Topology<sub_elem_dim> &topology)
{
  auto &ref_elem = phys_elem_.get_ref_space_element();
  ref_space_handler_.template init_cache<sub_elem_dim>(ref_elem);

  auto &map_elem = phys_elem_.get_map_element();
  phys_domain_handler_.template init_cache<sub_elem_dim>(map_elem);


  using RefSpHndlr = ReferenceElementHandler<dim_,range_,rank_>;
  const auto &grid_handler = dynamic_cast<RefSpHndlr &>(ref_space_handler_).get_grid_handler();

  auto &all_sub_elems_cache = phys_elem_.get_all_sub_elems_cache();
  if (all_sub_elems_cache == nullptr)
  {
    using VCache = typename PhysSpace::ElementAccessor::parent_t::Cache;

    using Cache = AllSubElementsCache<VCache>;
    all_sub_elems_cache = std::make_shared<Cache>();
  }

//  const auto n_basis = ref_elem.get_max_num_basis();//ref_elem.get_num_basis(DofProperties::active);
  const auto n_basis = ref_elem.get_basis_offset()[BaseSpace::n_components];

  const auto n_points = grid_handler.template get_num_points<sub_elem_dim>();
  for (auto &s_id: UnitElement<dim>::template elems_ids<sub_elem_dim>())
  {
    auto &sub_elem_cache = all_sub_elems_cache->template get_sub_elem_cache<sub_elem_dim>(s_id);
    sub_elem_cache.resize(flags_[sub_elem_dim], n_points, n_basis);
  }
}

template<int dim_,int range_,int rank_,int codim_>
void
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
init_cache(SpaceElement<dim_,codim_,range_,rank_> &sp_elem,
           const topology_variant &topology)
{
  using PhysElem = PhysicalSpaceElement<dim_,range_,rank_,codim_>;
  PhysElem *as_phys_elem = dynamic_cast<PhysElem *>(&sp_elem);
  Assert(as_phys_elem != nullptr,ExcNullPtr());

  auto init_cache_dispatcher =
    InitCacheDispatcher(flags_,*ref_space_handler_,phys_domain_handler_,*as_phys_elem);
  boost::apply_visitor(init_cache_dispatcher,topology);
}




template<int dim_,int range_,int rank_,int codim_>
template<int sub_elem_dim>
void
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
FillCacheDispatcher::
operator()(const Topology<sub_elem_dim> &topology)
{
  auto &ref_elem = phys_elem_.get_ref_space_element();
  ref_space_handler_.template fill_cache<sub_elem_dim>(ref_elem, sub_elem_id_);

  auto &map_elem = phys_elem_.get_map_element();
  phys_domain_handler_.template fill_cache<sub_elem_dim>(map_elem, sub_elem_id_);

  auto &all_sub_elems_cache = phys_elem_.get_all_sub_elems_cache();
  Assert(all_sub_elems_cache != nullptr, ExcNullPtr());
  auto &sub_elem_cache = all_sub_elems_cache->template get_sub_elem_cache<sub_elem_dim>(sub_elem_id_);


//    using std::cref;
  if (sub_elem_cache.template status_fill<_Value>())
  {
    PushFwd::template
    transform_0<RefBasis::range,RefBasis::rank,sub_elem_dim>(sub_elem_id_,ref_elem,map_elem,phys_elem_);

    sub_elem_cache.template set_status_filled<_Value>(true);
  }
  if (sub_elem_cache.template status_fill<_Gradient>())
  {
    PushFwd::template
    transform_1<RefBasis::range,RefBasis::rank,sub_elem_dim>(sub_elem_id_,ref_elem,map_elem,phys_elem_);

    sub_elem_cache.template set_status_filled<_Gradient>(true);
  }
  if (sub_elem_cache.template status_fill<_Hessian>())
  {
    PushFwd::template
    transform_2<RefBasis::range,RefBasis::rank,sub_elem_dim>(sub_elem_id_,ref_elem,map_elem,phys_elem_);

    sub_elem_cache.template set_status_filled<_Hessian>(true);
  }
  if (sub_elem_cache.template status_fill<_Divergence>())
  {
    eval_divergences_from_gradients(
      sub_elem_cache.template get_data<_Gradient>(),
      sub_elem_cache.template get_data<_Divergence>());

    sub_elem_cache.template set_status_filled<_Divergence>(true);
  }

  sub_elem_cache.set_filled(true);
}

template<int dim_,int range_,int rank_,int codim_>
void
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
fill_cache(SpaceElement<dim_,codim_,range_,rank_> &sp_elem,
           const topology_variant &topology,
           const int sub_elem_id)
{
  using PhysElem = PhysicalSpaceElement<dim_,range_,rank_,codim_>;
  PhysElem *as_phys_elem = dynamic_cast<PhysElem *>(&sp_elem);
  Assert(as_phys_elem != nullptr,ExcNullPtr());

  auto fill_cache_dispatcher =
    FillCacheDispatcher(sub_elem_id,*ref_space_handler_,phys_domain_handler_,*as_phys_elem);
  boost::apply_visitor(fill_cache_dispatcher,topology);
}
#endif




template<int dim_,int range_,int rank_,int codim_>
auto
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
print_info(LogStream &out) const -> void
{
  ref_space_handler_->print_info(out);
  //  PFCache::print_info(out);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/phys_space_element_handler.inst>
