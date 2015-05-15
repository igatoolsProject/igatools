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



template<int dim,int range,int rank,int codim>
PhysSpaceElementHandler<dim,range,rank,codim>::
PhysSpaceElementHandler(std::shared_ptr<const PhysSpace> space)
    :
    base_t(space),
    ref_space_handler_(space->get_reference_space()->get_elem_handler()),
    push_fwd_(space->get_map_func())
{}

template<int dim,int range,int rank,int codim>
auto
PhysSpaceElementHandler<dim,range,rank,codim>::
create(std::shared_ptr<const PhysSpace> space) -> std::shared_ptr<self_t>
{
    Assert(space != nullptr,ExcNullPtr());
    return std::shared_ptr<self_t>(new self_t(space));
}



template<int dim,int range,int rank,int codim>
template<int k>
void
PhysSpaceElementHandler<dim,range,rank,codim>::
reset(const ValueFlags flag, const Quadrature<k> &eval_pts)
{
    /*
    ref_space_handler_->reset(space_to_ref_flag(PhysSpace::PushForwardType::type, flag), eval_pts);

    ValueFlags map_flags;
    TransformationFlags transf_flags;
    space_to_pf_flag(flag,map_flags, transf_flags);

    push_fwd_.template reset<k>(map_flags, transf_flags, eval_pts);
    flags_[k] = flag;
    //*/

    const std::set<int> elems_id =
        this->get_space()->get_grid()->get_elements_id();

    this->reset_selected_elements(
        flag,
        eval_pts,
        SafeSTLVector<int>(elems_id.begin(),elems_id.end()));
}


template<int dim,int range,int rank,int codim>
template<int k>
void
PhysSpaceElementHandler<dim,range,rank,codim>::
reset_selected_elements(
    const ValueFlags &flag,
    const Quadrature<k> &eval_pts,
    const SafeSTLVector<Index> &elements_flat_id)
{
    ref_space_handler_->
    reset_selected_elements(space_to_ref_flag(PhysSpace::PushForwardType::type, flag), eval_pts, elements_flat_id);


    ValueFlags map_flags;
    TransformationFlags transf_flags;
    space_to_pf_flag(flag,map_flags, transf_flags);

    push_fwd_.template reset<k>(map_flags, transf_flags, eval_pts);
    flags_[k] = flag;
}


template<int dim,int range,int rank,int codim>
template<int k>
void
PhysSpaceElementHandler<dim,range,rank,codim>::
reset_one_element(
    const ValueFlags &flag,
    const Quadrature<k> &eval_pts,
    const int elem_flat_id)
{
    this->reset_selected_elements(flag,eval_pts,SafeSTLVector<int>(1,elem_flat_id));
}



template<int dim,int range,int rank,int codim>
template<int k>
void
PhysSpaceElementHandler<dim,range,rank,codim>::
init_cache(ElementAccessor &elem)
{
    auto &ref_elem = elem.get_ref_space_element();
    ref_space_handler_->template init_cache<k>(ref_elem);

    auto &push_fwd_elem = elem.get_push_forward_accessor();
    push_fwd_.template init_cache<k>(push_fwd_elem);

    auto &cache = elem.PhysSpace::ElementAccessor::parent_t::all_sub_elems_cache_;
    if (cache == nullptr)
    {
        using VCache = typename PhysSpace::ElementAccessor::parent_t::Cache;

        using Cache = LocalCache<VCache>;
        cache = std::make_shared<Cache>();
    }

    const auto n_basis = ref_elem.get_num_basis(DofProperties::active);
    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_sub_elem_cache<k>(s_id);
        const auto n_points = ref_space_handler_->template get_num_points<k>();

        s_cache.resize(flags_[k], n_points, n_basis);
    }
}



template<int dim,int range,int rank,int codim>
template<int k>
void
PhysSpaceElementHandler<dim,range,rank,codim>::
fill_cache(ElementAccessor &elem, const int j)
{
    auto &ref_elem = elem.get_ref_space_element();
    ref_space_handler_->template fill_cache<k>(ref_elem, j);

    auto &push_fwd_elem = elem.get_push_forward_accessor();
    push_fwd_.template fill_cache<k>(push_fwd_elem, j);

    auto &local_cache = elem.PhysSpace::ElementAccessor::parent_t::all_sub_elems_cache_;
    Assert(local_cache != nullptr, ExcNullPtr());
    auto &cache = local_cache->template get_sub_elem_cache<k>(j);

    using std::cref;
    if (cache.template status_fill<_Value>())
    {
        auto &result = cache.template get_data<_Value>();
        const auto &ref_values = ref_elem.template get_basis<_Value,k>(j,DofProperties::active);
        push_fwd_elem.template transform_0<RefSpace::range,RefSpace::rank>
        (ref_values, result);

        cache.template set_status_filled<_Value>(true);
    }
    if (cache.template status_fill<_Gradient>())
    {
        const auto &ref_values = ref_elem.template get_basis<   _Value,k>(j,DofProperties::active);
        const auto &ref_der_1  = ref_elem.template get_basis<_Gradient,k>(j,DofProperties::active);
        const auto &values = cache.template get_data<_Value>();
        push_fwd_elem.template transform_1<PhysSpace::range,PhysSpace::rank, k>(
            std::make_tuple(cref(ref_values), cref(ref_der_1)),
            values,
            cache.template get_data<_Gradient>(),
            j);

        cache.template set_status_filled<_Gradient>(true);
    }
    if (cache.template status_fill<_Hessian>())
    {
        const auto &ref_values = ref_elem.template get_basis<   _Value,k>(j,DofProperties::active);
        const auto &ref_der_1  = ref_elem.template get_basis<_Gradient,k>(j,DofProperties::active);
        const auto &ref_der_2  = ref_elem.template get_basis< _Hessian,k>(j,DofProperties::active);
        const auto &values = cache.template get_data<   _Value>();
        const auto &der_1  = cache.template get_data<_Gradient>();
        push_fwd_elem.template transform_2<PhysSpace::range,PhysSpace::rank, k>(
            std::make_tuple(cref(ref_values), cref(ref_der_1), cref(ref_der_2)),
            std::make_tuple(cref(values),cref(der_1)),
            cache.template get_data<_Hessian>(),
            j);

        cache.template set_status_filled<_Hessian>(true);
    }
    if (cache.template status_fill<_Divergence>())
    {
        eval_divergences_from_gradients(
            cache.template get_data<_Gradient>(),
            cache.template get_data<_Divergence>());

        cache.template set_status_filled<_Divergence>(true);
    }

    cache.set_filled(true);

}

#if 0
template<int dim, int codim,int range,int rank>
auto
PhysSpaceElementHandler<dim,range,rank,codim>::
fill_element_cache(ElementAccessor &elem) -> void
{
    auto &ref_elem = elem.get_ref_space_element();
    RefPhysSpaceElementHandler::fill_element_cache(ref_elem);
    PFCache::fill_element(elem);

    auto &cache = elem.PhysSpace::ElementAccessor::parent_t::all_sub_elems_cache_;
    auto &elem_cache = cache->template get_sub_elem_cache<0>(0);

    if (elem_cache.flags_handler_.fill_values())
    {
        const auto &ref_values = ref_elem.template get_basis_ders<0,0>(0);
        elem.template transform_0<RefSpace::range,RefSpace::rank>
        (ref_values, elem_cache.template get_data<0>());

        elem_cache.flags_handler_.set_values_filled(true);
    }


    if (elem_cache.flags_handler_.fill_gradients())
    {
        const auto &ref_values = ref_elem.template get_basis_ders<0,0>(0);
        const auto &ref_der_1  = ref_elem.template get_basis_ders<0,1>(0);
        const auto &values = elem_cache.template get_data<0>();
        elem.template transform_1<PhysSpace::range,PhysSpace::rank>
        (std::make_tuple(ref_values, ref_der_1), values,
        elem_cache.template get_data<1>());

        elem_cache.flags_handler_.set_gradients_filled(true);
    }

#if 0

    if (cache.flags_handler_.fill_hessians())
    {
        if (transformation_type == Transformation::h_grad)
        {
            ValueTable<typename RefElemAccessor::Value> dummy;
            PfElemAccessor::
            template transform_hessians<PhysSpace::range,PhysSpace::rank>(
                dummy,
                ref_space_element_accessor_.get_basis_gradients(topology_id),
                ref_space_element_accessor_.get_basis_hessians(topology_id),
                elem_cache.template get_datas<2>(),
                topology_id);

        }
        else
        {
            Assert(false,ExcNotImplemented());
            AssertThrow(false,ExcNotImplemented());

        }
        cache.flags_handler_.set_hessians_filled(true);
    }

    if (cache.flags_handler_.fill_divergences())
    {
        Assert(cache.flags_handler_.gradients_filled(),
               ExcMessage("Divergence requires gradient to be filled."));

        auto D1  = cache.D1phi_.begin();
        auto div = cache.div_phi_.begin();
        auto end = cache.D1phi_.end();
        for (; D1 != end; ++D1, ++div)
            *div = trace(*D1);

        cache.flags_handler_.set_divergences_filled(true);
        //Assert(false,ExcNotImplemented());
        //AssertThrow(false,ExcNotImplemented());
    }

#endif

    elem_cache.set_filled(true);
}
#endif





template<int dim,int range,int rank,int codim>
auto
PhysSpaceElementHandler<dim,range,rank,codim>::
print_info(LogStream &out) const -> void
{
    ref_space_handler_->print_info(out);
    //  PFCache::print_info(out);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/phys_space_element_handler.inst>
