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

#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/base/exceptions.h>

using std::array;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template< class PhysSpace >
PhysicalSpaceElement<PhysSpace>::
PhysicalSpaceElement(const std::shared_ptr<ContainerType> phys_space,
                     const Index index)
    :
    parent_t(phys_space,index),
//    PfElemAccessor(phys_space->get_map_func(), index),
    ref_space_element_accessor_(phys_space->get_reference_space()->create_element(index)),
    push_fwd_element_(shared_ptr<PfElemAccessor>(new PfElemAccessor(phys_space->get_map_func(), index)))
{
    Assert(ref_space_element_accessor_ != nullptr, ExcNullPtr());
    Assert(push_fwd_element_ != nullptr, ExcNullPtr());
}



template< class PhysSpace >
PhysicalSpaceElement<PhysSpace>::
PhysicalSpaceElement(const std::shared_ptr<ContainerType> phys_space,
                     const TensorIndex<dim> &index)
    :
    PhysicalSpaceElement(phys_space,phys_space->get_grid()->tensor_to_flat(index))
{}


template< class PhysSpace >
PhysicalSpaceElement<PhysSpace>::
PhysicalSpaceElement(const PhysicalSpaceElement<PhysSpace> &in,
                     const CopyPolicy &copy_policy)
    :
    parent_t(in,copy_policy)
    //,
//    PfElemAccessor(in,copy_policy)
{
    if (copy_policy == CopyPolicy::shallow)
    {
        ref_space_element_accessor_ = in.ref_space_element_accessor_;
        push_fwd_element_ = in.push_fwd_element_;
    }
    else
    {
        ref_space_element_accessor_ =
            shared_ptr<RefElemAccessor>(new RefElemAccessor(*in.ref_space_element_accessor_));
        push_fwd_element_ =
            shared_ptr<PfElemAccessor>(new PfElemAccessor(*in.push_fwd_element_));
    }

    Assert(false,ExcNotImplemented());
}




template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
copy_from(const PhysicalSpaceElement<PhysSpace> &element,
          const CopyPolicy &copy_policy)
{
    Assert(false,ExcNotImplemented());
//    SpaceElementAccessor<PhysSpace>::copy_from(element,copy_policy);
//
//    PhysSpace::PushForwardType::ElementAccessor::copy_from(element,copy_policy);
//
//    if (copy_policy == CopyPolicy::deep)
//        ref_space_element_accessor_->deep_copy_from(element.ref_space_element_accessor_);
//    else if (copy_policy == CopyPolicy::shallow)
//        ref_space_element_accessor_->deep_copy_from(element.ref_space_element_accessor_);
//    else
//    {
//        Assert(false,ExcNotImplemented());
//    }
}

template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
deep_copy_from(const PhysicalSpaceElement<PhysSpace> &element)
{
    Assert(false,ExcNotImplemented());
    //this->copy_from(element,CopyPolicy::deep);
}


template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
shallow_copy_from(const PhysicalSpaceElement<PhysSpace> &element)
{
    Assert(false,ExcNotImplemented());
//    this->copy_from(element,CopyPolicy::shallow);
}

#if 0
template< class PhysSpace >
ValueFlags
PhysicalSpaceElement<PhysSpace>::
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
PhysicalSpaceElement<PhysSpace>::
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

    if (contains(fill_flag , ValueFlags::measure))
        reference_space_accessor_fill_flags |= ValueFlags::measure;


    return reference_space_accessor_fill_flags;
}



template< class PhysSpace >
ValueFlags
PhysicalSpaceElement<PhysSpace>::
get_push_forward_accessor_fill_flags(const ValueFlags fill_flag) const
{
    const ValueFlags common_flag =
        ValueFlags::point|
        ValueFlags::map_value|
        ValueFlags::map_gradient|
        ValueFlags::map_hessian|
        ValueFlags::map_inv_gradient|
        ValueFlags::map_inv_hessian|
        ValueFlags::measure|
        ValueFlags::w_measure|
        ValueFlags::face_point|
        ValueFlags::map_face_value|
        ValueFlags::map_face_gradient|
        ValueFlags::map_face_hessian|
        ValueFlags::map_face_inv_gradient|
        ValueFlags::map_face_inv_hessian|
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
#endif

#if 0
template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
init_cache(const ValueFlags fill_flag,
           const QuadratureType &quad)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());

    const ValueFlags ref_sp_flag =
        get_reference_space_accessor_fill_flags(fill_flag);
    // TODO (pauletti, Sep 12, 2014): fix next line
    // ref_space_element_accessor_->init_cache(ref_sp_flag, quad);

    //const ValueFlags pf_flag = get_push_forward_accessor_fill_flags(fill_flag);
    //PfElemAccessor::init_cache(pf_flag, quad);


    // TODO (pauletti, Sep 12, 2014): fix next line
    // this->reset_element_and_faces_cache(fill_flag, quad);
}



template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
init_face_cache(const Index face_id,
                const ValueFlags fill_flag,
                const QuadratureFaceType &quad)
{
    AssertThrow(false,ExcNotImplemented());
}



template< class PhysSpace >
void PhysicalSpaceElement<PhysSpace>::
fill_cache(const TopologyId<dim> &topology_id)
{
    auto &cache = parent_t::get_values_cache(topology_id);

    Assert(cache.is_initialized(), ExcNotInitialized());

    //TODO: remove this if
    if (topology_id.is_element())
    {
        PfElemAccessor::fill_cache();
        // TODO (pauletti, Sep 12, 2014): fix next line
        // ref_space_element_accessor_->fill_cache();
    }
    else
    {
        //TODO: implement fill_cache in PushForwardElementAccessor
        // and RefSpaceElementAccessor accepting TopologyId
        PfElemAccessor::fill_face_cache(topology_id.get_id());
        // TODO (pauletti, Sep 12, 2014): fix next line
        //ref_space_element_accessor_->fill_face_cache(topology_id.get_id());
    }

    if (cache.flags_handler_.fill_values())
    {
        PfElemAccessor::
        template transform_values<RefSpace::range,RefSpace::rank>(
            ref_space_element_accessor_->get_basis_values(topology_id),
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
                ref_space_element_accessor_->get_basis_gradients(topology_id),
                cache.D1phi_,
                topology_id);
        }
        else
        {
            PfElemAccessor::
            template transform_gradients<PhysSpace::range,PhysSpace::rank>(
                ref_space_element_accessor_->get_basis_values(topology_id),
                ref_space_element_accessor_->get_basis_gradients(topology_id),
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
                ref_space_element_accessor_->get_basis_gradients(topology_id),
                ref_space_element_accessor_->get_basis_hessians(topology_id),
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

    cache.set_filled(true);
}



template< class PhysSpace >
void PhysicalSpaceElement<PhysSpace>::
fill_face_cache(const Index face_id)
{
    this->fill_cache(FaceTopology<dim>(face_id));
}

#endif

#if 0
template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_point(const Index qp,const TopologyId<dim> &topology_id) const -> const PhysPoint &
{
//    Assert(this->get_values_cache(topology_id).is_filled(), ExcCacheNotFilled());
    return (PfElemAccessor::get_map_values(topology_id))[qp];
}



template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
evaluate_field(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Value >
{
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(), local_coefs.size()));

    auto field_hat = ref_space_element_accessor_->evaluate_field(local_coefs,topology_id);

    ValueVector<Value> field(field_hat.size());

    PfElemAccessor::template
    transform_values<PhysSpace::range,PhysSpace::rank>(field_hat, field);

    return field;
}



template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
evaluate_field_gradients(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Derivative<1> >
{
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(), local_coefs.size()));

    auto D1field_hat = ref_space_element_accessor_->evaluate_field_gradients(local_coefs,topology_id);

    const auto n_quad_points = D1field_hat.size();

    ValueVector< typename RefElemAccessor::Value > D0field_hat(n_quad_points);
    if (transformation_type != Transformation::h_grad)
        D0field_hat = ref_space_element_accessor_->evaluate_field(local_coefs,topology_id);

    ValueVector< Derivative<1> > D1field(n_quad_points);
    PfElemAccessor::
    template transform_gradients<PhysSpace::range,PhysSpace::rank>(
        D0field_hat, D1field_hat, D1field);

    return D1field;
}



template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
evaluate_field_hessians(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Derivative<2> >
{
    AssertThrow(false,ExcNotImplemented());
    ValueVector< Derivative<2> > D2field;

    return D2field;
}









template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_points(const TopologyId<dim> &topology_id) const ->
const ValueVector< typename Mapping<dim, codim>::Value > &
{
//    Assert(this->get_values_cache(topology_id).is_filled(), ExcCacheNotFilled());
    return PfElemAccessor::get_map_values(topology_id);
}

template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_face_points(const Index face_id) const ->
const ValueVector< typename Mapping<dim, codim>::Value > &
{
    return this->get_points(FaceTopology<dim>(face_id));
}

/*
template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
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
PhysicalSpaceElement<PhysSpace>::
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
PhysicalSpaceElement<PhysSpace>::
is_boundary() const
{
    return PfElemAccessor::is_boundary();
}
//*/


template< class PhysSpace >
bool
PhysicalSpaceElement<PhysSpace>::
is_boundary(const Index face) const
{
    return PfElemAccessor::is_boundary(face);
}
//*/

#endif
template< class PhysSpace >
Index
PhysicalSpaceElement<PhysSpace>::
get_flat_index() const
{
    return parent_t::get_flat_index();
}

template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_tensor_index() const -> TensorIndex<dim>
{
    return parent_t::get_tensor_index();
}

#if 0
template< class PhysSpace >
bool
PhysicalSpaceElement<PhysSpace>::
jump(const TensorIndex<dim> &increment)
{
    const bool jump_grid_accessor = parent_t::jump(increment);

    const bool jump_push_fwd_accessor = PfElemAccessor::jump(increment);

    const bool jump_ref_space_accessor = ref_space_element_accessor_->jump(increment);

    return jump_grid_accessor && jump_push_fwd_accessor && jump_ref_space_accessor;
}
#endif

template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
move_to(const Index flat_index)
{
    this->as_cartesian_grid_element_accessor().move_to(flat_index);
    ref_space_element_accessor_->move_to(flat_index);
    push_fwd_element_->move_to(flat_index);
}

/*
template< class PhysSpace >
bool
PhysicalSpaceElement<PhysSpace>::
operator==(const PhysicalSpaceElement <PhysSpace> &a) const
{
    return this->as_cartesian_grid_element_accessor() == a.as_cartesian_grid_element_accessor();
}



template< class PhysSpace >
bool
PhysicalSpaceElement<PhysSpace>::
operator!=(const PhysicalSpaceElement <PhysSpace> &a) const
{
    return this->as_cartesian_grid_element_accessor() != a.as_cartesian_grid_element_accessor();
}

template< class PhysSpace >
bool
PhysicalSpaceElement<PhysSpace>::
operator>(const PhysicalSpaceElement <PhysSpace> &a) const
{
    return this->as_cartesian_grid_element_accessor() > a.as_cartesian_grid_element_accessor();
}

template< class PhysSpace >
bool
PhysicalSpaceElement<PhysSpace>::
operator<(const PhysicalSpaceElement <PhysSpace> &a) const
{
    return this->as_cartesian_grid_element_accessor() < a.as_cartesian_grid_element_accessor();
}
//*/

template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_ref_space_accessor() const -> const RefElemAccessor &
{
    return *ref_space_element_accessor_;
}

template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_ref_space_accessor() -> RefElemAccessor &
{
    return *ref_space_element_accessor_;
}

template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_grid() const -> const std::shared_ptr<const CartesianGrid<dim> >
{
    return this->get_ref_space_accessor().get_grid();
}

template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_push_forward_accessor() const -> const PfElemAccessor &
{
    return *push_fwd_element_;
}

template< class PhysSpace >
auto
PhysicalSpaceElement<PhysSpace>::
get_push_forward_accessor() -> PfElemAccessor &
{
    return *push_fwd_element_;
}


#if 0
template< class PhysSpace >
template <int deriv_order>
auto
PhysicalSpaceElement<PhysSpace>::
evaluate_basis_derivatives_at_points(const ValueVector<RefPoint> &points) const ->
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
        phi_hat = ref_space_element_accessor_->evaluate_basis_values_at_points(points);

    ValueTable<ref_gradients_t> D1phi_hat;
    if (contains(ref_space_flags,ValueFlags::gradient))
        D1phi_hat = ref_space_element_accessor_->evaluate_basis_gradients_at_points(points);

    ValueTable<ref_hessians_t> D2phi_hat;
    if (contains(ref_space_flags,ValueFlags::hessian))
        D2phi_hat = ref_space_element_accessor_->evaluate_basis_hessians_at_points(points);
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

#endif

template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
print_info(LogStream &out) const
{
    out.begin_item("Reference space:");
    ref_space_element_accessor_->print_info(out);
    out.end_item();

    out.begin_item("Pushforward:");
    push_fwd_element_->print_info(out);
    out.end_item();
}

template< class PhysSpace >
void
PhysicalSpaceElement<PhysSpace>::
print_cache_info(LogStream &out) const
{
    out.begin_item("Reference space:");
    ref_space_element_accessor_->print_cache_info(out);
    out.end_item();

    out.begin_item("Pushforward:");
    push_fwd_element_->print_cache_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

//#include <igatools/basis_functions/physical_space_element.inst>
