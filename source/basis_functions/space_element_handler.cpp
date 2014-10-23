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
NewValueFlags
get_reference_space_accessor_fill_flags(const NewValueFlags fill_flag,
                                        const Transformation transformation_type)
{
    bool fill_values = false;
    bool fill_gradients = false;
    bool fill_hessians = false;
//    bool fill_face_values = false;
//    bool fill_face_gradients = false;
//    bool fill_face_hessians = false;

    if (contains(fill_flag , NewValueFlags::value))
        fill_values = true;

    if (contains(fill_flag , NewValueFlags::gradient))
        fill_gradients = true;

    if (contains(fill_flag , NewValueFlags::hessian))
        fill_hessians = true;

//    if (contains(fill_flag , NewValueFlags::face_value))
//        fill_face_values = true;
//
//    if (contains(fill_flag , NewValueFlags::face_gradient))
//        fill_face_gradients = true;
//
//    if (contains(fill_flag , NewValueFlags::face_hessian))
//        fill_face_hessians = true;


    bool fill_D0_phi_hat = false;
    bool fill_D1_phi_hat = false;
    bool fill_D2_phi_hat = false;
//    bool fill_face_D0_phi_hat = false;
//    bool fill_face_D1_phi_hat = false;
//    bool fill_face_D2_phi_hat = false;
    if (transformation_type == Transformation::h_grad)
    {
        fill_D0_phi_hat = fill_values;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
//        fill_face_D0_phi_hat = fill_face_values;
//        fill_face_D1_phi_hat = fill_face_gradients || fill_face_hessians;
//        fill_face_D2_phi_hat = fill_face_hessians;
    }
    else if (transformation_type == Transformation::h_div  ||
             transformation_type == Transformation::h_curl ||
             transformation_type == Transformation::l_2)
    {
        fill_D0_phi_hat = fill_values || fill_gradients || fill_hessians;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
//        fill_face_D0_phi_hat = fill_face_values || fill_face_gradients || fill_face_hessians;
//        fill_face_D1_phi_hat = fill_face_gradients || fill_face_hessians;
//        fill_face_D2_phi_hat = fill_face_hessians;
    }


    NewValueFlags reference_space_accessor_fill_flags = NewValueFlags::none;
    if (fill_D0_phi_hat)
        reference_space_accessor_fill_flags |= NewValueFlags::value;

    if (fill_D1_phi_hat)
        reference_space_accessor_fill_flags |= NewValueFlags::gradient;

    if (fill_D2_phi_hat)
        reference_space_accessor_fill_flags |= NewValueFlags::hessian;

//    if (fill_face_D0_phi_hat)
//        reference_space_accessor_fill_flags |= NewValueFlags::face_value;
//
//    if (fill_face_D1_phi_hat)
//        reference_space_accessor_fill_flags |= NewValueFlags::face_gradient;
//
//    if (fill_face_D2_phi_hat)
//        reference_space_accessor_fill_flags |= NewValueFlags::face_hessian;

    if (contains(fill_flag , NewValueFlags::measure))
        reference_space_accessor_fill_flags |= NewValueFlags::measure;


    return reference_space_accessor_fill_flags;
}


NewValueFlags
get_push_forward_accessor_fill_flags(const NewValueFlags fill_flag)
{
    const NewValueFlags common_flag =
        NewValueFlags::point|
        NewValueFlags::value|
        NewValueFlags::gradient|
        NewValueFlags::hessian|
        NewValueFlags::measure|
        NewValueFlags::w_measure;

//        NewValueFlags::face_point|
//        NewValueFlags::map_face_value|
//        NewValueFlags::map_face_gradient|
//        NewValueFlags::map_face_hessian|
//        NewValueFlags::face_w_measure|
//        NewValueFlags::face_normal;

    NewValueFlags pf_flags = fill_flag & common_flag;

    if (contains(fill_flag , NewValueFlags::value))
        pf_flags |= NewValueFlags::tran_value;

    if (contains(fill_flag , NewValueFlags::gradient))
        pf_flags |= NewValueFlags::tran_gradient;

    if (contains(fill_flag , NewValueFlags::hessian))
        pf_flags |= NewValueFlags::tran_hessian;

    return pf_flags;
}

};



template<class PhysSpace>
SpaceElementHandler<PhysSpace>::
SpaceElementHandler(std::shared_ptr<const PhysSpace> space,
                    const NewValueFlags flag,
                    const Quadrature<dim> &quad)
    :
    RefSpaceElementHandler(space->get_reference_space(), get_reference_space_accessor_fill_flags(flag, PhysSpace::PushForwardType::type), quad),
    PFCache(space->get_map_func(), get_push_forward_accessor_fill_flags(flag), quad),
    space_(space),
    flags_ {flag, flag},
       quad_(quad)
{}



template<class PhysSpace>
auto
SpaceElementHandler<PhysSpace>::
init_element_cache(ElementAccessor &elem) -> void
{
    RefSpaceElementHandler::init_element_cache(elem.get_ref_space_accessor());
    PFCache::init_element(elem);

    auto &cache = elem.PhysSpace::ElementAccessor::parent_t::local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename PhysSpace::ElementAccessor::parent_t::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    auto n_basis = space_->get_num_all_element_basis();
    auto &elem_cache = cache->template get_value_cache<0>(0);
    elem_cache.resize(std::get<0>(flags_), quad_, n_basis);

    for (auto &f : RefSpaceElementHandler::faces)
    {
        auto &face_cache = cache->template get_value_cache<1>(f);
        face_cache.resize(std::get<1>(flags_), quad_.collapse_to_face(f), n_basis);
    }
}



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



template<class PhysSpace>
auto
SpaceElementHandler<PhysSpace>::
init_element_cache(ElementIterator &elem) -> void
{
    init_element_cache(elem.get_accessor());
}


template<class PhysSpace>
auto
SpaceElementHandler<PhysSpace>::
fill_element_cache(ElementIterator &elem) -> void
{
    fill_element_cache(elem.get_accessor());
}



template<class PhysSpace>
auto
SpaceElementHandler<PhysSpace>::
print_info(LogStream &out) const -> void
{
    RefSpaceElementHandler::print_info(out);
    //  PFCache::print_info(out);
}


IGA_NAMESPACE_CLOSE

//#include <igatools/basis_functions/space_element_handler.inst>
