//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#include <igatools/basis_functions/space_element_handler.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace
{
auto
space_to_ref_flag(const Transformation type, const NewValueFlags flags)
-> NewValueFlags
{
    NewValueFlags ref_flag = NewValueFlags::none;

    bool fill_values = false;
    bool fill_gradients = false;
    bool fill_hessians = false;

    if (contains(flags , NewValueFlags::value))
        fill_values = true;
    if (contains(flags , NewValueFlags::gradient))
        fill_gradients = true;
    if (contains(flags , NewValueFlags::hessian))
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
        ref_flag |= NewValueFlags::value;
    if (fill_D1_phi_hat)
        ref_flag |= NewValueFlags::gradient;
    if (fill_D2_phi_hat)
        ref_flag |= NewValueFlags::hessian;

    return ref_flag;
}


auto
space_to_pf_flag(const NewValueFlags flags)
{
    NewValueFlags transfer_flag =
        NewValueFlags::measure |
        NewValueFlags::w_measure |
        NewValueFlags::outer_normal |
        NewValueFlags::boundary_normal |
        NewValueFlags::point;

    NewValueFlags pf_flag = flags & transfer_flag;

    if (contains(flags , NewValueFlags::value))
        pf_flag |= NewValueFlags::tran_value;

    if (contains(flags , NewValueFlags::gradient))
        pf_flag |= NewValueFlags::tran_gradient;

    if (contains(flags , NewValueFlags::hessian))
        pf_flag |= NewValueFlags::tran_hessian;

    return pf_flag;
}

};



template<class PhysSpace>
SpaceElementHandler<PhysSpace>::
SpaceElementHandler(std::shared_ptr<const PhysSpace> space)
    :
    RefSpaceElementHandler(space->get_reference_space()),
    PFCache(space->get_map_func()),
    space_(space)
{}



template<class PhysSpace>
template<int k>
void
SpaceElementHandler<PhysSpace>::
reset(const NewValueFlags flag, const Quadrature<k> &quad)
{
    RefSpaceElementHandler::template reset<k>(space_to_ref_flag(type, flag), quad);
    PFCache::template reset<k>(space_to_pf_flag(flag), quad);
    flags_[k] = flag;
}



template<class PhysSpace>
template<int k>
void
SpaceElementHandler<PhysSpace>::
init_cache(ElementAccessor &elem)
{
    auto &ref_elem = elem.get_ref_space_accessor();
    RefSpaceElementHandler::template init_cache<k>(ref_elem);
    PFCache::template init_cache<k>(elem);

    auto &cache = elem.PhysSpace::ElementAccessor::parent_t::local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename PhysSpace::ElementAccessor::parent_t::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    //auto n_basis = space_->get_num_all_element_basis();
    const auto n_basis = ref_elem.get_num_basis();

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_value_cache<k>(s_id);
        const auto n_points = this->template get_num_points<k>();

        s_cache.resize(flags_[k], n_points, n_basis);
    }
}



template<class PhysSpace>
template<int k>
void
SpaceElementHandler<PhysSpace>::
fill_cache(ElementAccessor &elem, const int j)
{
    auto &ref_elem = elem.get_ref_space_accessor();
    RefSpaceElementHandler::template fill_cache<k>(ref_elem, j);
    PFCache::template fill_cache<k>(elem, j);

    auto &local_cache = elem.PhysSpace::ElementAccessor::parent_t::local_cache_;
    Assert(local_cache != nullptr, ExcNullPtr());
    auto &cache =  local_cache->template get_value_cache<k>(j);

    auto &flags = cache.flags_handler_;

    if (flags.fill_values())
    {
        auto &result = cache.template get_der<0>();
        const auto &ref_values = ref_elem.template get_values<0,k>(j);
        elem.template transform_0<RefSpace::range,RefSpace::rank>
        (ref_values, result);

        flags.set_values_filled(true);
    }
    if (flags.fill_gradients())
    {
        const auto &ref_values = ref_elem.template get_values<0,k>(j);
        const auto &ref_der_1  = ref_elem.template get_values<1,k>(j);
        const auto &values = cache.template get_der<0>();
        elem.template transform_1<PhysSpace::range,PhysSpace::rank, k>
        (std::make_tuple(ref_values, ref_der_1), values,
         cache.template get_der<1>(), j);

        flags.set_gradients_filled(true);
    }
    if (flags.fill_hessians())
    {
        const auto &ref_values = ref_elem.template get_values<0,k>(j);
        const auto &ref_der_1  = ref_elem.template get_values<1,k>(j);
        const auto &ref_der_2  = ref_elem.template get_values<2,k>(j);
        const auto &values = cache.template get_der<0>();
        const auto &der_1  = cache.template get_der<1>();
        elem.template transform_2<PhysSpace::range,PhysSpace::rank, k>
        (std::make_tuple(ref_values, ref_der_1, ref_der_2),
         std::make_tuple(values,der_1),
         cache.template get_der<2>(), j);

        flags.set_hessians_filled(true);
    }

    cache.set_filled(true);

}

#if 0
template<class PhysSpace>
auto
SpaceElementHandler<PhysSpace>::
fill_element_cache(ElementAccessor &elem) -> void
{
    auto &ref_elem = elem.get_ref_space_accessor();
    RefSpaceElementHandler::fill_element_cache(ref_elem);
    PFCache::fill_element(elem);

    auto &cache = elem.PhysSpace::ElementAccessor::parent_t::local_cache_;
    auto &elem_cache = cache->template get_value_cache<0>(0);

    if (elem_cache.flags_handler_.fill_values())
    {
        const auto &ref_values = ref_elem.template get_basis_ders<0,0>(0);
        elem.template transform_0<RefSpace::range,RefSpace::rank>
        (ref_values, elem_cache.template get_der<0>());

        elem_cache.flags_handler_.set_values_filled(true);
    }


    if (elem_cache.flags_handler_.fill_gradients())
    {
        const auto &ref_values = ref_elem.template get_basis_ders<0,0>(0);
        const auto &ref_der_1  = ref_elem.template get_basis_ders<0,1>(0);
        const auto &values = elem_cache.template get_der<0>();
        elem.template transform_1<PhysSpace::range,PhysSpace::rank>
        (std::make_tuple(ref_values, ref_der_1), values,
        elem_cache.template get_der<1>());

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
                elem_cache.template get_ders<2>(),
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





template<class PhysSpace>
auto
SpaceElementHandler<PhysSpace>::
print_info(LogStream &out) const -> void
{
    RefSpaceElementHandler::print_info(out);
    //  PFCache::print_info(out);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_element_handler.inst>
