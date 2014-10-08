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

#include <igatools/basis_functions/space_uniform_quad_cache.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace
{
ValueFlags
get_reference_space_accessor_fill_flags(const ValueFlags fill_flag,
                                        const Transformation transformation_type)
{
    bool fill_values = false;
    bool fill_gradients = false;
    bool fill_hessians = false;
    bool fill_face_values = false;
    bool fill_face_gradients = false;
    bool fill_face_hessians = false;

    if (contains(fill_flag , ValueFlags::value))
        fill_values = true;

    if (contains(fill_flag , ValueFlags::gradient))
        fill_gradients = true;

    if (contains(fill_flag , ValueFlags::hessian))
        fill_hessians = true;

    if (contains(fill_flag , ValueFlags::face_value))
        fill_face_values = true;

    if (contains(fill_flag , ValueFlags::face_gradient))
        fill_face_gradients = true;

    if (contains(fill_flag , ValueFlags::face_hessian))
        fill_face_hessians = true;


    bool fill_D0_phi_hat = false;
    bool fill_D1_phi_hat = false;
    bool fill_D2_phi_hat = false;
    bool fill_face_D0_phi_hat = false;
    bool fill_face_D1_phi_hat = false;
    bool fill_face_D2_phi_hat = false;
    if (transformation_type == Transformation::h_grad)
    {
        fill_D0_phi_hat = fill_values;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
        fill_face_D0_phi_hat = fill_face_values;
        fill_face_D1_phi_hat = fill_face_gradients || fill_face_hessians;
        fill_face_D2_phi_hat = fill_face_hessians;
    }
    else if (transformation_type == Transformation::h_div  ||
             transformation_type == Transformation::h_curl ||
             transformation_type == Transformation::l_2)
    {
        fill_D0_phi_hat = fill_values || fill_gradients || fill_hessians;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
        fill_face_D0_phi_hat = fill_face_values || fill_face_gradients || fill_face_hessians;
        fill_face_D1_phi_hat = fill_face_gradients || fill_face_hessians;
        fill_face_D2_phi_hat = fill_face_hessians;
    }


    ValueFlags reference_space_accessor_fill_flags = ValueFlags::none;
    if (fill_D0_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::value;

    if (fill_D1_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::gradient;

    if (fill_D2_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::hessian;

    if (fill_face_D0_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::face_value;

    if (fill_face_D1_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::face_gradient;

    if (fill_face_D2_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::face_hessian;

    if (contains(fill_flag , ValueFlags::measure))
        reference_space_accessor_fill_flags |= ValueFlags::measure;


    return reference_space_accessor_fill_flags;
}


ValueFlags
get_push_forward_accessor_fill_flags(const ValueFlags fill_flag)
{
    const ValueFlags common_flag =
        ValueFlags::point|
        ValueFlags::map_value|
        ValueFlags::map_gradient|
        ValueFlags::map_hessian|
        ValueFlags::measure|
        ValueFlags::w_measure|
        ValueFlags::face_point|
        ValueFlags::map_face_value|
        ValueFlags::map_face_gradient|
        ValueFlags::map_face_hessian|
        ValueFlags::face_w_measure|
        ValueFlags::face_normal;

    ValueFlags pf_flags = fill_flag & common_flag;

    if (contains(fill_flag , ValueFlags::value) || contains(fill_flag , ValueFlags::face_value))
        pf_flags |= ValueFlags::tran_value;

    if (contains(fill_flag , ValueFlags::gradient) || contains(fill_flag , ValueFlags::face_gradient))
        pf_flags |= ValueFlags::tran_gradient;

    if (contains(fill_flag , ValueFlags::hessian) || contains(fill_flag , ValueFlags::face_hessian))
        pf_flags |= ValueFlags::tran_hessian;

    return pf_flags;
}

};



template<class PhysSpace>
SpaceUniformQuadCache<PhysSpace>::
SpaceUniformQuadCache(std::shared_ptr<const PhysSpace> space,
                          const ValueFlags flag,
                          const Quadrature<dim> &quad)
        :
        RefSpaceCache(space->get_reference_space(), get_reference_space_accessor_fill_flags(flag, PhysSpace::PushForwardType::transformation_type), quad),
        PFCache(space->get_push_forward(), get_push_forward_accessor_fill_flags(flag), quad),
        space_(space),
        flags_(flag),
        quad_(quad)
{}



template<class PhysSpace>
auto
SpaceUniformQuadCache<PhysSpace>::
init_element_cache(ElementAccessor &elem) -> void
{
    RefSpaceCache::init_element_cache(elem.get_ref_space_accessor());
    PFCache::init_element_cache(elem);
    auto &cache = elem.PhysSpace::ElementAccessor::parent_t::local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename PhysSpace::ElementAccessor::parent_t::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    auto n_basis = space_->get_num_all_element_basis();
    auto &elem_cache = cache->elem_values_;
    elem_cache.resize(flags_, quad_, n_basis);

    //        auto &face_cache = cache->face_values_;
    //        for (auto f: base_t::faces)
    //            face_cache[f].resize(face_flags_, quad_, n_basis, f);

}



template<class PhysSpace>
auto
SpaceUniformQuadCache<PhysSpace>::
fill_element_cache(ElementAccessor &elem) -> void
{
    auto &ref_elem = elem.get_ref_space_accessor();
    RefSpaceCache::fill_element_cache(ref_elem);
    PFCache::fill_element_cache(elem);

    auto &cache = elem.get_elem_cache();

    if (cache.flags_handler_.fill_values())
    {
        elem.template transform_values<RefSpace::range,RefSpace::rank>
        (ref_elem.get_basis_values(), cache.phi_);

        cache.flags_handler_.set_values_filled(true);
    }

    cache.set_filled(true);
}



template<class PhysSpace>
auto
SpaceUniformQuadCache<PhysSpace>::
init_element_cache(ElementIterator &elem) -> void
{
    init_element_cache(elem.get_accessor());
}


template<class PhysSpace>
auto
SpaceUniformQuadCache<PhysSpace>::
fill_element_cache(ElementIterator &elem) -> void
{
    fill_element_cache(elem.get_accessor());
}



template<class PhysSpace>
auto
SpaceUniformQuadCache<PhysSpace>::
print_info(LogStream &out) const -> void
{
    RefSpaceCache::print_info(out);
    PFCache::print_info(out);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_uniform_quad_cache.inst>
