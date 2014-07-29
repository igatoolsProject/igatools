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

#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/base/exceptions.h>

using std::array;
using std::vector;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template< class PhysSpace >
PhysicalSpaceElementAccessor<PhysSpace>::
PhysicalSpaceElementAccessor(const std::shared_ptr<ContainerType> phys_space,
                             const Index index)
    :
    SpaceElementAccessor<
    PhysicalSpaceElementAccessor<PhysSpace>,
    PhysSpace,
    PhysSpace::RefSpace::dim,
    PhysSpace::PushForwardType::codim,
    PhysSpace::RefSpace::range,
    PhysSpace::RefSpace::rank>(phys_space,index),
    PfElemAccessor(phys_space->get_push_forward(), index),
    ref_space_element_accessor_(phys_space->get_reference_space(),index)
{}







template< class PhysSpace >
ValueFlags
PhysicalSpaceElementAccessor<PhysSpace>::
get_face_flags(const ValueFlags fill_flag) const
{
    ValueFlags face_fill_flag = ValueFlags::none ;

    if (contains(fill_flag , ValueFlags::face_value))
        face_fill_flag |= ValueFlags::value ;

    if (contains(fill_flag , ValueFlags::face_gradient))
        face_fill_flag |= ValueFlags::gradient ;

    if (contains(fill_flag , ValueFlags::face_hessian))
        face_fill_flag |= ValueFlags::hessian ;

    return face_fill_flag ;
}



template< class PhysSpace >
ValueFlags
PhysicalSpaceElementAccessor<PhysSpace>::
get_reference_space_accessor_fill_flags(const ValueFlags fill_flag) const
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


    return reference_space_accessor_fill_flags;
}



template< class PhysSpace >
ValueFlags
PhysicalSpaceElementAccessor<PhysSpace>::
get_push_forward_accessor_fill_flags(const ValueFlags fill_flag) const
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



template< class PhysSpace >
void
PhysicalSpaceElementAccessor<PhysSpace>::
init_values(const ValueFlags fill_flag,
            const QuadratureType &quad)
{
    const ValueFlags ref_sp_flag =
        get_reference_space_accessor_fill_flags(fill_flag);
    ref_space_element_accessor_.init_values(ref_sp_flag, quad);

    const ValueFlags pf_flag = get_push_forward_accessor_fill_flags(fill_flag);
    PfElemAccessor::init_values(pf_flag, quad);


    this->reset_element_and_faces_cache(fill_flag, quad);
}



template< class PhysSpace >
void
PhysicalSpaceElementAccessor<PhysSpace>::
init_face_values(const Index face_id,
                 const ValueFlags fill_flag,
                 const QuadratureFaceType &quad)
{
    AssertThrow(false,ExcNotImplemented());
}



template< class PhysSpace >
void PhysicalSpaceElementAccessor<PhysSpace>::
fill_values(const TopologyId<dim> &topology_id)
{
    auto &cache = parent_t::get_values_cache(topology_id);

    Assert(cache.is_initialized(), ExcNotInitialized());

    //TODO: remove this if
    if (topology_id.is_element())
    {
        PfElemAccessor::fill_values();
        ref_space_element_accessor_.fill_values();
    }
    else
    {
        //TODO: implement fill_values in PushForwardElementAccessor
        // and RefSpaceElementAccessor accepting TopologyId
        PfElemAccessor::fill_face_values(topology_id.get_id());
        ref_space_element_accessor_.fill_face_values(topology_id.get_id());
    }

    if (cache.flags_handler_.fill_values())
    {
        PfElemAccessor::
        template transform_values<RefSpace::range,RefSpace::rank>(
            ref_space_element_accessor_.get_basis_values(topology_id),
            cache.phi_,
            topology_id);

        cache.flags_handler_.set_values_filled(true);
    }

    if (cache.flags_handler_.fill_gradients())
    {
        if (transformation_type == Transformation::h_grad)
        {
            ValueTable<typename RefElemAccessor::Value> dummy;
            PfElemAccessor::
            template transform_gradients<PhysSpace::range,PhysSpace::rank>(
                dummy,
                ref_space_element_accessor_.get_basis_gradients(topology_id),
                cache.D1phi_,
                topology_id);
        }
        else
        {
            PfElemAccessor::
            template transform_gradients<PhysSpace::range,PhysSpace::rank>(
                ref_space_element_accessor_.get_basis_values(topology_id),
                ref_space_element_accessor_.get_basis_gradients(topology_id),
                cache.D1phi_,
                topology_id);
        }
        cache.flags_handler_.set_gradients_filled(true);
    }

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
                cache.D2phi_,
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
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
    }

    cache.set_filled(true);
}



template< class PhysSpace >
void PhysicalSpaceElementAccessor<PhysSpace>::
fill_face_values(const Index face_id)
{
    this->fill_values(FaceTopology<dim>(face_id));
}


template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_point(const Index qp,const TopologyId<dim> &topology_id) const -> const PhysPoint &
{
//    Assert(this->get_values_cache(topology_id).is_filled(), ExcCacheNotFilled());
    return (PfElemAccessor::get_map_values(topology_id))[qp];
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
evaluate_field(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Value >
{
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(), local_coefs.size()));

    auto field_hat = ref_space_element_accessor_.evaluate_field(local_coefs,topology_id);

    ValueVector<Value> field(field_hat.size());

    PfElemAccessor::template
    transform_values<PhysSpace::range,PhysSpace::rank>(field_hat, field);

    return field;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
evaluate_field_gradients(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Derivative<1> >
{
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(), local_coefs.size()));

    auto D1field_hat = ref_space_element_accessor_.evaluate_field_gradients(local_coefs,topology_id);

    const auto n_quad_points = D1field_hat.size();

    ValueVector< typename RefElemAccessor::Value > D0field_hat(n_quad_points);
    if (transformation_type != Transformation::h_grad)
        D0field_hat = ref_space_element_accessor_.evaluate_field(local_coefs,topology_id);

    ValueVector< Derivative<1> > D1field(n_quad_points);
    PfElemAccessor::
    template transform_gradients<PhysSpace::range,PhysSpace::rank>(
        D0field_hat, D1field_hat, D1field);

    return D1field;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
evaluate_field_hessians(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Derivative<2> >
{
    AssertThrow(false,ExcNotImplemented());
    ValueVector< Derivative<2> > D2field;

    return D2field;
}









template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_points(const TopologyId<dim> &topology_id) const ->
const ValueVector< typename Mapping<dim, codim>::Value > &
{
//    Assert(this->get_values_cache(topology_id).is_filled(), ExcCacheNotFilled());
    return PfElemAccessor::get_map_values(topology_id);
}

template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_face_points(const Index face_id) const ->
const ValueVector< typename Mapping<dim, codim>::Value > &
{
    return this->get_points(FaceTopology<dim>(face_id));
}

/*
template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_map_gradients(const TopologyId<dim> &topology_id) const ->
const ValueVector< typename Mapping<dim, codim>::Gradient > &
{
//    Assert(this->get_values_cache(topology_id).is_filled(), ExcCacheNotFilled());
    return PfElemAccessor::get_gradients(topology_id);
}
//*/


/*
template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_face_normals(const Index face_id) const ->
const ValueVector< typename Mapping<dim, codim>::Value > &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
//    Assert(this->face_values_[face_id].is_filled(), ExcCacheNotFilled());
    return PfElemAccessor::get_face_normals(face_id);
}
//*/

template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
is_boundary() const
{
    return PfElemAccessor::is_boundary();
}
//*/


template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
is_boundary(const Index face) const
{
    return PfElemAccessor::is_boundary(face);
}
//*/

template< class PhysSpace >
Index
PhysicalSpaceElementAccessor<PhysSpace>::
get_flat_index() const
{
    return parent_t::get_flat_index();
}
//*/

template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_physical_space() const -> std::shared_ptr<const PhysSpace>
{
    return this->space_;
}



template< class PhysSpace >
void
PhysicalSpaceElementAccessor<PhysSpace>::
operator++()
{
    parent_t::operator++();
    PfElemAccessor::operator++();
    ++ref_space_element_accessor_;
}



template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
operator==(const PhysicalSpaceElementAccessor <PhysSpace> &a) const
{
    return this->as_cartesian_grid_element_accessor() == a.as_cartesian_grid_element_accessor();
}



template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
operator!=(const PhysicalSpaceElementAccessor <PhysSpace> &a) const
{
    return this->as_cartesian_grid_element_accessor() != a.as_cartesian_grid_element_accessor();
}

template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_ref_space_accessor() const -> const RefElemAccessor &
{
    return ref_space_element_accessor_;
}

template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_push_forward_accessor() const -> const PfElemAccessor &
{
    return static_cast<const PfElemAccessor &>(*this);
}




template< class PhysSpace >
template <int deriv_order>
auto
PhysicalSpaceElementAccessor<PhysSpace>::
evaluate_basis_derivatives_at_points(const std::vector<RefPoint> &points) const ->
ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
{
    Assert(deriv_order >= 0 && deriv_order <= 2,ExcIndexRange(deriv_order,0,2));

    const Size n_basis  = this->get_num_basis();
    const Size n_points = points.size();


    ValueFlags phys_space_flags;
    if (deriv_order == 0)
    {
        phys_space_flags = ValueFlags::value;
    }
    else if (deriv_order == 1)
    {
        phys_space_flags = ValueFlags::gradient;
    }
    else if (deriv_order == 2)
    {
        phys_space_flags = ValueFlags::hessian;
    }


    //---------------------------------------------------------------------------------------------
    // evaluation of the basis function values (or derivatives) using the reference space --- begin
    const ValueFlags ref_space_flags = get_reference_space_accessor_fill_flags(phys_space_flags);

    using ref_values_t = typename RefElemAccessor::Value;
    using ref_gradients_t = typename RefElemAccessor::template Derivative<1>;
    using ref_hessians_t = typename RefElemAccessor::template Derivative<2>;

    ValueTable<ref_values_t> phi_hat;
    if (contains(ref_space_flags,ValueFlags::value))
        phi_hat = ref_space_element_accessor_.evaluate_basis_values_at_points(points);

    ValueTable<ref_gradients_t> D1phi_hat;
    if (contains(ref_space_flags,ValueFlags::gradient))
        D1phi_hat = ref_space_element_accessor_.evaluate_basis_gradients_at_points(points);

    ValueTable<ref_hessians_t> D2phi_hat;
    if (contains(ref_space_flags,ValueFlags::hessian))
        D2phi_hat = ref_space_element_accessor_.evaluate_basis_hessians_at_points(points);
    // evaluation of the basis function values (or derivatives) using the reference space --- end
    //---------------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------------------------
    // basis function push forwarding from ref. space to phys. space --- begin
    ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > > transformed_basis(n_basis,n_points);
    PfElemAccessor::
    template transform_basis_derivatives_at_points<PhysSpace::range,PhysSpace::rank>(
        points,
        phi_hat,
        D1phi_hat,
        D2phi_hat,
        transformed_basis);
    // basis function push forwarding from ref. space to phys. space --- end
    //---------------------------------------------------------------------------------------------


    return transformed_basis;
}



template< class PhysSpace >
void
PhysicalSpaceElementAccessor<PhysSpace>::
print_info(LogStream &out, const VerbosityLevel verbosity_level) const
{
    using std::endl ;

    std::string tab = "   ";

    out << "PhysicalSpaceElementAccessor info:" << endl;
    out.push(tab);

    ref_space_element_accessor_.print_info(out,verbosity_level);
    PfElemAccessor::print_info(out,verbosity_level);

    out.pop();
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_space_element_accessor.inst>
