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
auto
space_to_ref_flag(const Transformation type, const ValueFlags flags)
-> ValueFlags
{
    ValueFlags ref_flag = ValueFlags::none;

    bool fill_values = false;
    bool fill_gradients = false;
    bool fill_hessians = false;

    if (contains(flags , ValueFlags::value))
        fill_values = true;
    if (contains(flags , ValueFlags::gradient))
        fill_gradients = true;
    if (contains(flags , ValueFlags::hessian))
        fill_hessians = true;

    bool fill_D0_phi_hat = false;
    bool fill_D1_phi_hat = false;
    bool fill_D2_phi_hat = false;

    if (type == Transformation::h_grad)
    {
        fill_D0_phi_hat = fill_values;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
    }
    else if (type == Transformation::h_div  ||
    type == Transformation::h_curl ||
    type == Transformation::l_2)
    {
        fill_D0_phi_hat = fill_values || fill_gradients || fill_hessians;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
    }

    if (fill_D0_phi_hat)
        ref_flag |= ValueFlags::value;
    if (fill_D1_phi_hat)
        ref_flag |= ValueFlags::gradient;
    if (fill_D2_phi_hat)
        ref_flag |= ValueFlags::hessian;

    return ref_flag;
}


void
space_to_pf_flag(const ValueFlags flags, ValueFlags &map_flags, TransformationFlags &transf_flags)
{
    ValueFlags transfer_flag =
        ValueFlags::measure |
        ValueFlags::w_measure |
        ValueFlags::outer_normal |
        ValueFlags::boundary_normal |
        ValueFlags::point;

    map_flags = ValueFlags::none;
    map_flags = flags & transfer_flag;



    transf_flags = TransformationFlags::tran_none;
    if (contains(flags , ValueFlags::value))
        transf_flags |= TransformationFlags::tran_value;

    if (contains(flags , ValueFlags::gradient))
        transf_flags |= TransformationFlags::tran_gradient;

    if (contains(flags , ValueFlags::hessian))
        transf_flags |= TransformationFlags::tran_hessian;
}

};



template<int dim_,int range_,int rank_,int codim_>
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
PhysSpaceElementHandler(std::shared_ptr<const PhysSpace> space)
    :
    base_t(space),
    ref_space_handler_(space->get_reference_space()->get_elem_handler()),
    push_fwd_(space->get_map_func())
{}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
create(std::shared_ptr<const PhysSpace> space) -> std::shared_ptr<self_t>
{
    Assert(space != nullptr,ExcNullPtr());
    return std::shared_ptr<self_t>(new self_t(space));
}




template<int dim_,int range_,int rank_,int codim_>
template<int sub_elem_dim>
void
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
ResetDispatcher::
operator()(const Quadrature<sub_elem_dim> &quad)
{
    flags_[sub_elem_dim] = flag_in_;
    ref_space_handler_.reset_selected_elements(
        space_to_ref_flag(PhysSpace::PushForwardType::type, flag_in_),
        quad,
        elements_flat_id_);


    ValueFlags map_flags;
    TransformationFlags transf_flags;
    space_to_pf_flag(flag_in_,map_flags, transf_flags);

    push_fwd_.template reset<sub_elem_dim>(map_flags,transf_flags,quad);
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
        ResetDispatcher(flag,elements_flat_id,*ref_space_handler_,push_fwd_,flags_);
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

    auto &push_fwd_elem = phys_elem_.get_push_forward_accessor();
    push_fwd_.template init_cache<sub_elem_dim>(push_fwd_elem);


    using RefSpHndlr = ReferenceElementHandler<dim_,range_,rank_>;
    const auto &grid_handler = dynamic_cast<RefSpHndlr &>(ref_space_handler_).get_grid_handler();

    auto &all_sub_elems_cache = phys_elem_.get_all_sub_elems_cache();
    if (all_sub_elems_cache == nullptr)
    {
        using VCache = typename PhysSpace::ElementAccessor::parent_t::Cache;

        using Cache = AllSubElementsCache<VCache>;
        all_sub_elems_cache = std::make_shared<Cache>();
    }

    const auto n_basis = ref_elem.get_max_num_basis();//ref_elem.get_num_basis(DofProperties::active);
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
        InitCacheDispatcher(flags_,*ref_space_handler_,push_fwd_,*as_phys_elem);
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

    auto &push_fwd_elem = phys_elem_.get_push_forward_accessor();
    push_fwd_.template fill_cache<sub_elem_dim>(push_fwd_elem, sub_elem_id_);

    auto &all_sub_elems_cache = phys_elem_.get_all_sub_elems_cache();
    Assert(all_sub_elems_cache != nullptr, ExcNullPtr());
    auto &sub_elem_cache = all_sub_elems_cache->template get_sub_elem_cache<sub_elem_dim>(sub_elem_id_);

    using std::cref;
    if (sub_elem_cache.template status_fill<_Value>())
    {
        auto &result = sub_elem_cache.template get_data<_Value>();
        const auto &ref_values = ref_elem.template get_basis<_Value,sub_elem_dim>(sub_elem_id_,DofProperties::active);
        push_fwd_elem.template transform_0<RefSpace::range,RefSpace::rank>
        (ref_values, result);

        sub_elem_cache.template set_status_filled<_Value>(true);
    }
    if (sub_elem_cache.template status_fill<_Gradient>())
    {
        const auto &ref_values = ref_elem.template get_basis<   _Value,sub_elem_dim>(sub_elem_id_,DofProperties::active);
        const auto &ref_der_1  = ref_elem.template get_basis<_Gradient,sub_elem_dim>(sub_elem_id_,DofProperties::active);
        const auto &values = sub_elem_cache.template get_data<_Value>();
        push_fwd_elem.template transform_1<PhysSpace::range,PhysSpace::rank,sub_elem_dim>(
            std::make_tuple(cref(ref_values), cref(ref_der_1)),
            values,
            sub_elem_cache.template get_data<_Gradient>(),
            sub_elem_id_);

        sub_elem_cache.template set_status_filled<_Gradient>(true);
    }
    if (sub_elem_cache.template status_fill<_Hessian>())
    {
        const auto &ref_values = ref_elem.template get_basis<   _Value,sub_elem_dim>(sub_elem_id_,DofProperties::active);
        const auto &ref_der_1  = ref_elem.template get_basis<_Gradient,sub_elem_dim>(sub_elem_id_,DofProperties::active);
        const auto &ref_der_2  = ref_elem.template get_basis< _Hessian,sub_elem_dim>(sub_elem_id_,DofProperties::active);
        const auto &values = sub_elem_cache.template get_data<   _Value>();
        const auto &der_1  = sub_elem_cache.template get_data<_Gradient>();
        push_fwd_elem.template transform_2<PhysSpace::range,PhysSpace::rank, sub_elem_dim>(
            std::make_tuple(cref(ref_values), cref(ref_der_1), cref(ref_der_2)),
            std::make_tuple(cref(values),cref(der_1)),
            sub_elem_cache.template get_data<_Hessian>(),
            sub_elem_id_);

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
        FillCacheDispatcher(sub_elem_id,*ref_space_handler_,push_fwd_,*as_phys_elem);
    boost::apply_visitor(fill_cache_dispatcher,topology);
}





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
